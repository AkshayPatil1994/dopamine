/* SPDX-FileCopyrightText: 2026 DOPAMINE contributors */
/* SPDX-License-Identifier: AGPL-3.0-only */

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Runtime-selected pressure preconditioner: boomeramg (default) or gamg. */
static char g_pc_solver_type[64] = "boomeramg";

#ifdef FVM_HAVE_PETSC
#include <petscksp.h>

typedef struct {
  int ready;
  int n_cells_local;
  int n_cells_global;
  int n_faces;
  int ksp_tolerances_set;
  int use_direct_solution_copy;
  /* global_n_dir is derived from n_dirichlet_rows, which is fixed by mesh/BC topology (set once
   * in pressure_solver_init) and never varies across solves for a given handle; cache it instead
   * of an MPI_Allreduce + null-space rebuild on every call. nullsp_attached avoids rebuilding the
   * null space when the attach/detach decision hasn't changed since the last solve. */
  int global_n_dir;
  int global_n_dir_valid;
  int nullsp_attached;  /* -1 = unknown, 0 = detached, 1 = attached */
  PetscInt ref_row;
  /* Local offset of the reference row in the owned block; -1 if not local. */
  PetscInt ref_local_offset;
  PetscInt ownership_start;
  PetscInt last_max_iterations;
  PetscReal last_rtol;
  PetscReal last_atol;
  Mat matrix;
  Vec rhs_vector;
  Vec solution_vector;
  KSP ksp;
  IS sol_from_is;
  Vec sol_local;
  VecScatter sol_scatter;
  PetscInt *sol_rows;
  PetscInt *sol_local_offsets;
  PetscScalar *rhs_values;
  PetscScalar *sol_values;
} fvm_petsc_pressure_ctx_t;

static fvm_petsc_pressure_ctx_t g_pressure_ctx = {
  0, 0, 0, 0, 0, 0,
  0, 0, -1,
  0, -1, 0, 0, 0.0, 0.0,
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL
};

static int fvm_petsc_gpu_backend_enabled(void) {
#if defined(PETSC_HAVE_CUDA)
  const char *env = getenv("DOPAMINE_PETSC_USE_GPU");
  return !(env != NULL && env[0] == '0' && env[1] == '\0');
#else
  return 0;
#endif
}

static int fvm_petsc_gpu_aware_mpi_enabled(void) {
#if defined(PETSC_HAVE_CUDA)
  const char *env = getenv("DOPAMINE_PETSC_USE_GPU_AWARE_MPI");
  return !(env != NULL && env[0] == '0' && env[1] == '\0');
#else
  return 0;
#endif
}

static PetscErrorCode fvm_petsc_apply_gpu_defaults(void) {
#if defined(PETSC_HAVE_CUDA)
  if (!fvm_petsc_gpu_backend_enabled()) return PETSC_SUCCESS;

  /* Mat/Vec CUDA types are set directly via MatSetType(MATAIJCUSPARSE)/VecSetType(VECCUDA)
   * right after creation (below), not through MatSetFromOptions/VecSetFromOptions -- so
   * -mat_type/-vec_type options database entries here would never be consumed. PETSc warns
   * about exactly that ("unused options") at exit if set; leave type selection to the
   * direct API calls only. -use_gpu_aware_mpi is set at PetscInitialize time instead, in
   * fvm_petsc_backend_initialize -- see the comment there. */

  if (strncmp(g_pc_solver_type, "boomeramg", 9) == 0) {
    PetscCall(PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_strong_threshold", "0.25"));
    /* hypre's default coarsening/interpolation selection can pick host-only algorithms
     * (e.g. FF1 interpolation) even when the matrix is device-resident, segfaulting in
     * host code (hypre_BoomerAMGBuildFF1Interp) because the CSR arrays it walks live on
     * the GPU. Pin to the coarsening/interpolation/relaxation choices hypre's CUDA path
     * actually supports. */
    PetscCall(PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_coarsen_type", "PMIS"));
    PetscCall(PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_interp_type", "ext+i"));
    PetscCall(PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_relax_type_all", "l1scaled-Jacobi"));
  }
#endif
  return PETSC_SUCCESS;
}

static void fvm_petsc_pressure_ctx_reset(fvm_petsc_pressure_ctx_t *ctx) {
  ctx->ready = 0;
  ctx->n_cells_local = 0;
  ctx->n_cells_global = 0;
  ctx->n_faces = 0;
  ctx->ksp_tolerances_set = 0;
  ctx->use_direct_solution_copy = 0;
  ctx->global_n_dir = 0;
  ctx->global_n_dir_valid = 0;
  ctx->nullsp_attached = -1;
  ctx->ref_row = 0;
  ctx->ref_local_offset = -1;
  ctx->ownership_start = 0;
  ctx->last_max_iterations = 0;
  ctx->last_rtol = 0.0;
  ctx->last_atol = 0.0;
  ctx->matrix = NULL;
  ctx->rhs_vector = NULL;
  ctx->solution_vector = NULL;
  ctx->ksp = NULL;
  ctx->sol_from_is = NULL;
  ctx->sol_local = NULL;
  ctx->sol_scatter = NULL;
  ctx->sol_rows = NULL;
  ctx->sol_local_offsets = NULL;
  ctx->rhs_values = NULL;
  ctx->sol_values = NULL;
}

static int fvm_petsc_pressure_ctx_destroy(fvm_petsc_pressure_ctx_t *ctx) {
  PetscErrorCode ierr;

  if (ctx->sol_scatter != NULL) {
    ierr = VecScatterDestroy(&ctx->sol_scatter);
    if (ierr != PETSC_SUCCESS) return -1;
  }
  if (ctx->sol_local != NULL) {
    ierr = VecDestroy(&ctx->sol_local);
    if (ierr != PETSC_SUCCESS) return -2;
  }
  if (ctx->sol_from_is != NULL) {
    ierr = ISDestroy(&ctx->sol_from_is);
    if (ierr != PETSC_SUCCESS) return -3;
  }
  if (ctx->sol_rows != NULL) {
    ierr = PetscFree(ctx->sol_rows);
    if (ierr != PETSC_SUCCESS) return -4;
    ctx->sol_rows = NULL;
  }
  if (ctx->sol_local_offsets != NULL) {
    ierr = PetscFree(ctx->sol_local_offsets);
    if (ierr != PETSC_SUCCESS) return -4;
    ctx->sol_local_offsets = NULL;
  }
  if (ctx->rhs_values != NULL) {
    ierr = PetscFree(ctx->rhs_values);
    if (ierr != PETSC_SUCCESS) return -4;
    ctx->rhs_values = NULL;
  }
  if (ctx->sol_values != NULL) {
    ierr = PetscFree(ctx->sol_values);
    if (ierr != PETSC_SUCCESS) return -4;
    ctx->sol_values = NULL;
  }
  if (ctx->ksp != NULL) {
    ierr = KSPDestroy(&ctx->ksp);
    if (ierr != PETSC_SUCCESS) return -5;
  }
  if (ctx->solution_vector != NULL) {
    ierr = VecDestroy(&ctx->solution_vector);
    if (ierr != PETSC_SUCCESS) return -6;
  }
  if (ctx->rhs_vector != NULL) {
    ierr = VecDestroy(&ctx->rhs_vector);
    if (ierr != PETSC_SUCCESS) return -7;
  }
  if (ctx->matrix != NULL) {
    ierr = MatDestroy(&ctx->matrix);
    if (ierr != PETSC_SUCCESS) return -8;
  }

  fvm_petsc_pressure_ctx_reset(ctx);
  return 0;
}

static int fvm_petsc_pressure_ctx_build(fvm_petsc_pressure_ctx_t *ctx,
                                        const int n_cells_local,
                                        const int n_cells_global,
                                        const int n_faces,
                                        const int *local_to_global,
                                        const int *face_owner,
                                        const int *face_neighbor_local,
                                        const int *face_neighbor_global,
                                        const int *skip_face,
                                        const double *diagonal,
                                        const double *face_coeff,
                                        const int reference_cell,
                                        const int max_iterations,
                                        const double rtol,
                                        const double atol) {
  PetscErrorCode ierr;
  PetscInt ref_row;
  PetscInt row;
  PetscInt ksp_max_iterations;
  int f;
  PetscLogDouble t_setup_start = 0.0, t_setup_end = 0.0;

  PetscInt *d_nnz_arr = NULL;
  PetscInt *o_nnz_arr = NULL;
  Vec zero_vec = NULL;
  PetscInt ownership_start = 0;
  PetscInt ownership_end = 0;


  // --- Begin timing and user message at the true start of setup ---
  PetscPrintf(PETSC_COMM_WORLD,
              " Pressure solver: building pressure solver for %d cells "
              "(one-time setup) ...\n",
              n_cells_global);
  PetscTime(&t_setup_start);

  ref_row = (PetscInt) reference_cell - 1;
  if (ref_row < 0 || ref_row >= (PetscInt) n_cells_global) ref_row = 0;

  d_nnz_arr = (PetscInt *) calloc((size_t) n_cells_local, sizeof(PetscInt));
  o_nnz_arr = (PetscInt *) calloc((size_t) n_cells_local, sizeof(PetscInt));
  if (!d_nnz_arr || !o_nnz_arr) {
    free(d_nnz_arr);
    free(o_nnz_arr);
    return -10;
  }

  for (f = 0; f < n_cells_local; ++f) d_nnz_arr[f] = 1;
  for (f = 0; f < n_faces; ++f) {
    int owner_local_idx;
    int nbr_local_idx;
    if (skip_face[f] != 0) continue;
    if (face_neighbor_global[f] <= 0) continue;

    owner_local_idx = face_owner[f] - 1;
    if (owner_local_idx < 0 || owner_local_idx >= n_cells_local) continue;

    if (face_neighbor_local[f] > 0) {
      nbr_local_idx = face_neighbor_local[f] - 1;
      if (nbr_local_idx < 0 || nbr_local_idx >= n_cells_local) continue;
      d_nnz_arr[owner_local_idx] += 1;
      d_nnz_arr[nbr_local_idx] += 1;
    } else {
      o_nnz_arr[owner_local_idx] += 1;
    }
  }

  for (f = 0; f < n_cells_local; ++f)
    o_nnz_arr[f] = (o_nnz_arr[f] * 120) / 100;

  ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                      (PetscInt) n_cells_local, (PetscInt) n_cells_local,
                      (PetscInt) n_cells_global, (PetscInt) n_cells_global,
                      0, d_nnz_arr, 0, o_nnz_arr, &ctx->matrix);
  free(d_nnz_arr);
  free(o_nnz_arr);
  if (ierr != PETSC_SUCCESS) return -10;

#if defined(PETSC_HAVE_CUDA)
  if (fvm_petsc_gpu_backend_enabled()) {
    ierr = MatSetType(ctx->matrix, MATAIJCUSPARSE);
    if (ierr != PETSC_SUCCESS) return -10;
  }
#endif

  ierr = MatSetOption(ctx->matrix, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
  if (ierr != PETSC_SUCCESS) return -13;
  ierr = MatSetOption(ctx->matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  if (ierr != PETSC_SUCCESS) return -13;

  ierr = VecCreateMPI(PETSC_COMM_WORLD, (PetscInt) n_cells_local, (PetscInt) n_cells_global, &ctx->rhs_vector);
  if (ierr != PETSC_SUCCESS) return -11;
  ierr = VecDuplicate(ctx->rhs_vector, &ctx->solution_vector);
  if (ierr != PETSC_SUCCESS) return -12;

#if defined(PETSC_HAVE_CUDA)
  if (fvm_petsc_gpu_backend_enabled()) {
    ierr = VecSetType(ctx->rhs_vector, VECCUDA);
    if (ierr != PETSC_SUCCESS) return -11;
    ierr = VecSetType(ctx->solution_vector, VECCUDA);
    if (ierr != PETSC_SUCCESS) return -12;
  }
#endif

  for (row = 0; row < (PetscInt) n_cells_local; ++row) {
    PetscInt global_row = (PetscInt) local_to_global[row] - 1;
    PetscScalar diag_value = (PetscScalar) diagonal[row];
    if (!isfinite((double) diag_value) || fabs((double) diag_value) < 1.0e-20) diag_value = 1.0;
    ierr = MatSetValue(ctx->matrix, global_row, global_row, diag_value, ADD_VALUES);
    if (ierr != PETSC_SUCCESS) return -14;
  }

  for (row = 0; row < (PetscInt) n_faces; ++row) {
    PetscInt owner_local;
    PetscInt owner_global;
    PetscInt neighbor_global;
    PetscScalar coeff;

    if (skip_face[row] != 0) continue;
    if (face_neighbor_global[row] <= 0) continue;

    owner_local = (PetscInt) face_owner[row] - 1;
    if (owner_local < 0 || owner_local >= (PetscInt) n_cells_local) continue;
    owner_global = (PetscInt) local_to_global[owner_local] - 1;
    neighbor_global = (PetscInt) face_neighbor_global[row] - 1;
    coeff = (PetscScalar) face_coeff[row];
    if (!isfinite((double) coeff) || coeff == 0.0) continue;

    ierr = MatSetValue(ctx->matrix, owner_global, neighbor_global, -coeff, ADD_VALUES);
    if (ierr != PETSC_SUCCESS) return -16;
    if (face_neighbor_local[row] > 0) {
      ierr = MatSetValue(ctx->matrix, neighbor_global, owner_global, -coeff, ADD_VALUES);
      if (ierr != PETSC_SUCCESS) return -17;
    }
  }

  ierr = MatAssemblyBegin(ctx->matrix, MAT_FINAL_ASSEMBLY);
  if (ierr != PETSC_SUCCESS) return -18;
  ierr = MatAssemblyEnd(ctx->matrix, MAT_FINAL_ASSEMBLY);
  if (ierr != PETSC_SUCCESS) return -19;

  /* NOTE: Reference row zeroing is deferred to solve phase; null space or Dirichlet BCs handle gauge. */

  /* SPD hints: GAMG and BoomerAMG exploit these for smoother selection. */
  ierr = MatSetOption(ctx->matrix, MAT_SYMMETRIC, PETSC_TRUE);
  if (ierr != PETSC_SUCCESS) return -25;
  ierr = MatSetOption(ctx->matrix, MAT_SPD, PETSC_TRUE);
  if (ierr != PETSC_SUCCESS) return -25;

  ierr = KSPCreate(PETSC_COMM_WORLD, &ctx->ksp);
  if (ierr != PETSC_SUCCESS) return -26;
  ierr = fvm_petsc_apply_gpu_defaults();
  if (ierr != PETSC_SUCCESS) return -26;
  ierr = KSPSetOperators(ctx->ksp, ctx->matrix, ctx->matrix);
  if (ierr != PETSC_SUCCESS) return -27;
  ierr = KSPSetType(ctx->ksp, KSPCG);
  if (ierr != PETSC_SUCCESS) return -28;
  ierr = KSPSetNormType(ctx->ksp, KSP_NORM_UNPRECONDITIONED);
  if (ierr != PETSC_SUCCESS) return -28;
  {
    PC pc = NULL;
    ierr = KSPGetPC(ctx->ksp, &pc);
    if (ierr != PETSC_SUCCESS) return -29;

    if (strncmp(g_pc_solver_type, "boomeramg", 9) == 0) {
      /* BoomerAMG defaults; users can override via PETSc options. */
      ierr = PCSetType(pc, PCHYPRE);
      if (ierr != PETSC_SUCCESS) return -30;
      ierr = PCHYPRESetType(pc, "boomeramg");
      if (ierr != PETSC_SUCCESS) return -30;

      /* HMIS + FF1 are robust CPU defaults for FV pressure operators, but FF1 is a
       * host-only interpolation routine that segfaults if applied to a device-resident
       * matrix (hypre_BoomerAMGBuildFF1Interp walking GPU pointers as host memory).
       * fvm_petsc_apply_gpu_defaults() already set GPU-safe equivalents above; don't
       * clobber them here. */
      if (!fvm_petsc_gpu_backend_enabled()) {
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_coarsen_type",     "HMIS");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_interp_type",      "FF1");
        /* Standard strong-threshold default for 3-D pressure solves. */
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_strong_threshold", "0.5");
        /* Symmetric smoother across levels for the SPD pressure operator. */
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_relax_type_all",   "symmetric-sor/jacobi");
      }
      /* Keep aggressive coarsening off by default. */
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_agg_nl",           "0");
      /* Cap interpolation stencil width to control memory growth. */
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_P_max",            "4");
    } else {
      /* PETSc GAMG path for smaller cases or non-hypre builds. */
      PetscReal gamg_threshold[1] = {0.02};

      ierr = PCSetType(pc, PCGAMG);
      if (ierr != PETSC_SUCCESS) return -30;

      /* Prune weak couplings before coarsening. */
      ierr = PCGAMGSetThreshold(pc, gamg_threshold, 1);
      if (ierr != PETSC_SUCCESS) return -30;

      /* Stop coarsening near 200 equations/rank on coarse levels. */
      ierr = PCGAMGSetCoarseEqLim(pc, 200);
      if (ierr != PETSC_SUCCESS) return -30;

      /* Use one setup smoothing pass per level. */
      ierr = PCGAMGSetNSmooths(pc, 1);
      if (ierr != PETSC_SUCCESS) return -30;

      /* Repartitioning requires ParMETIS/Chaco/PTScotch; disable so GAMG
       * works with the default partitioner (no external graph library). */
      ierr = PCGAMGSetRepartition(pc, PETSC_FALSE);
      if (ierr != PETSC_SUCCESS) return -30;
    }
  }

  ksp_max_iterations = (PetscInt) ((max_iterations > 2000) ? max_iterations : 2000);
  ierr = KSPSetTolerances(ctx->ksp, (PetscReal) rtol, (PetscReal) atol, PETSC_DEFAULT, ksp_max_iterations);
  if (ierr != PETSC_SUCCESS) return -31;
  ctx->last_rtol = (PetscReal) rtol;
  ctx->last_atol = (PetscReal) atol;
  ctx->last_max_iterations = ksp_max_iterations;
  ctx->ksp_tolerances_set = 1;
  ierr = KSPSetInitialGuessNonzero(ctx->ksp, PETSC_TRUE);
  if (ierr != PETSC_SUCCESS) return -31;
  ierr = KSPSetFromOptions(ctx->ksp);
  if (ierr != PETSC_SUCCESS) return -32;

  /* Register constant null space so CG projects out {1,...,1} at every iteration. */
  {
    MatNullSpace nullsp = NULL;
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullsp);
    if (ierr == PETSC_SUCCESS) {
      ierr = MatSetNullSpace(ctx->matrix, nullsp);
      if (ierr != PETSC_SUCCESS)
        PetscPrintf(PETSC_COMM_WORLD,
          " WARNING: MatSetNullSpace failed (ierr=%d); ref-cell divU may be nonzero\n", ierr);
      MatNullSpaceDestroy(&nullsp);
    }
  }

  ierr = KSPSetUp(ctx->ksp);
  if (ierr != PETSC_SUCCESS) return -33;
  PetscTime(&t_setup_end);
  PetscPrintf(PETSC_COMM_WORLD,
              " Pressure solver: Setup complete (%.1f s)\n",
              (double)(t_setup_end - t_setup_start));
  if (fvm_petsc_gpu_backend_enabled()) {
    PetscPrintf(PETSC_COMM_WORLD,
                " Pressure solver: PETSc CUDA objects enabled%s\n",
                fvm_petsc_gpu_aware_mpi_enabled() ? " with GPU-aware MPI" : "");
  }

  ierr = PetscMalloc1((PetscInt) n_cells_local, &ctx->sol_rows);
  if (ierr != PETSC_SUCCESS) return -37;
  ierr = PetscMalloc1((PetscInt) n_cells_local, &ctx->rhs_values);
  if (ierr != PETSC_SUCCESS) return -41;
  ierr = PetscMalloc1((PetscInt) n_cells_local, &ctx->sol_values);
  if (ierr != PETSC_SUCCESS) return -42;
  for (row = 0; row < (PetscInt) n_cells_local; ++row) {
    PetscInt gid = (PetscInt) local_to_global[row];
    if (gid < 1 || gid > (PetscInt) n_cells_global) return -45;
    ctx->sol_rows[row] = gid - 1;
  }

  ierr = VecGetOwnershipRange(ctx->solution_vector, &ownership_start, &ownership_end);
  if (ierr != PETSC_SUCCESS) return -46;
  ctx->ownership_start = ownership_start;

  /* Check if all local rows are owned locally to enable the fast copy path. */
  ctx->use_direct_solution_copy = 1;
  ctx->ref_local_offset = -1;
  for (row = 0; row < (PetscInt) n_cells_local; ++row) {
    if (ctx->sol_rows[row] < ownership_start || ctx->sol_rows[row] >= ownership_end) {
      ctx->use_direct_solution_copy = 0;
      break;
    }
    if (ctx->sol_rows[row] == ctx->ref_row)
      ctx->ref_local_offset = ctx->sol_rows[row] - ownership_start;
  }

  if (!ctx->use_direct_solution_copy) {
    ierr = ISCreateGeneral(PETSC_COMM_SELF, (PetscInt) n_cells_local, ctx->sol_rows,
                           PETSC_COPY_VALUES, &ctx->sol_from_is);
    if (ierr != PETSC_SUCCESS) return -38;
    ierr = VecCreateSeq(PETSC_COMM_SELF, (PetscInt) n_cells_local, &ctx->sol_local);
    if (ierr != PETSC_SUCCESS) return -39;
    ierr = VecScatterCreate(ctx->solution_vector, ctx->sol_from_is, ctx->sol_local, NULL, &ctx->sol_scatter);
    if (ierr != PETSC_SUCCESS) return -40;
  }

  ctx->n_cells_local = n_cells_local;
  ctx->n_cells_global = n_cells_global;
  ctx->n_faces = n_faces;
  ctx->ref_row = ref_row;
  ctx->ready = 1;
  return 0;
}
#endif

/* -- Public functions available regardless of PETSc build -------------- */

// Set the preconditioner type for subsequent pressure solves.
 // Recognised names: "boomeramg" (default), "gamg".
void fvm_petsc_set_pressure_solver(const char *solver_name) {
  if (solver_name) {
    strncpy(g_pc_solver_type, solver_name, sizeof(g_pc_solver_type) - 1);
    g_pc_solver_type[sizeof(g_pc_solver_type) - 1] = '\0';
  }
}

// Insert one option string into PETSc's global options database.
 // Call this after fvm_petsc_backend_initialize so PETSc is initialised.
void fvm_petsc_insert_option(const char *opt_str) {
#ifdef FVM_HAVE_PETSC
  if (opt_str && opt_str[0] != '\0') {
    PetscOptionsInsertString(NULL, opt_str);
  }
#else
  (void) opt_str;
#endif
}

int fvm_petsc_backend_available(void) {
#ifdef FVM_HAVE_PETSC
  return 1;
#else
  return 0;
#endif
}

int fvm_petsc_backend_initialize(void) {
#ifdef FVM_HAVE_PETSC
  PetscBool is_initialized = PETSC_FALSE;
  char program_name[] = "dopamine";
  char no_x_option[] = "-nox";
  char no_signal_handler[] = "-no_signal_handler";
  char gpu_aware_mpi_option[] = "-use_gpu_aware_mpi";
  char gpu_aware_mpi_value_on[] = "1";
  char gpu_aware_mpi_value_off[] = "0";
  char *argv_values[6] = {program_name, no_x_option, no_signal_handler, NULL, NULL, NULL};
  char **argv = argv_values;
#if defined(PETSC_HAVE_CUDA)
  /* -use_gpu_aware_mpi is read during PETSc's device init, which runs as part of
   * PetscInitialize -- setting it afterwards via PetscOptionsSetValue (as
   * fvm_petsc_apply_gpu_defaults used to) is too late and PETSc flags it "unused" at exit.
   * Always pass it explicitly (1 or 0) rather than omitting it when off: PETSc's own
   * built-in default when built with --with-cuda-aware-mpi=1 is to assume GPU-aware MPI,
   * so simply not passing the flag silently ignores DOPAMINE_PETSC_USE_GPU_AWARE_MPI=0
   * and aborts anyway. Only pass it at all when PETSc actually has CUDA support compiled
   * in -- the option doesn't exist otherwise and PETSc would flag it "unused" at exit. */
  int use_gpu_aware_mpi = fvm_petsc_gpu_backend_enabled() && fvm_petsc_gpu_aware_mpi_enabled();
  int argc = 5;
  argv_values[3] = gpu_aware_mpi_option;
  argv_values[4] = use_gpu_aware_mpi ? gpu_aware_mpi_value_on : gpu_aware_mpi_value_off;
#else
  int argc = 3;
#endif
  if (PetscInitialized(&is_initialized) != PETSC_SUCCESS) {
    return -2;
  }
  if (!is_initialized) {
    if (PetscInitialize(&argc, &argv, NULL, NULL) != PETSC_SUCCESS) {
      return -3;
    }
  }
  return 0;
#else
  return -1;
#endif
}

int fvm_petsc_backend_finalize(void) {
#ifdef FVM_HAVE_PETSC
  PetscBool is_finalized = PETSC_FALSE;
  PetscBool is_initialized = PETSC_FALSE;
  int ctx_status;
  ctx_status = fvm_petsc_pressure_ctx_destroy(&g_pressure_ctx);
  if (ctx_status != 0) {
    return -20 + ctx_status;
  }
  if (PetscInitialized(&is_initialized) != PETSC_SUCCESS) {
    return -2;
  }
  if (PetscFinalized(&is_finalized) != PETSC_SUCCESS) {
    return -3;
  }
  if (is_initialized && !is_finalized) {
    if (PetscFinalize() != PETSC_SUCCESS) {
      return -4;
    }
  }
  return 0;
#else
  return -1;
#endif
}

int fvm_petsc_solve_pressure(const int n_cells_local,
                            const int n_cells_global,
                            const int n_faces,
                            const int *local_to_global,
                            const int *face_owner,
                            const int *face_neighbor_local,
                            const int *face_neighbor_global,
                            const int *skip_face,
                            const double *diagonal,
                            const double *rhs,
                            const double *face_coeff,
                            const int reference_cell,
                            const int max_iterations,
                            const double rtol,
                            const double atol,
                            const int *dirichlet_rows,
                            const double *dirichlet_values,
                            const int n_dirichlet_rows,
                            double *solution,
                            int *iteration_count,
                            double *residual_norm) {
#ifdef FVM_HAVE_PETSC
  PetscErrorCode ierr;
  PetscInt row;
  PetscInt ref_row_expected;
  PetscInt ksp_max_iterations;
  PetscReal current_rtol;
  PetscReal current_atol;
  PetscReal petsc_residual = 0.0;
  KSPConvergedReason reason;
  int status = 0;
  int ctx_status;

  ref_row_expected = (PetscInt) reference_cell - 1;
  if (ref_row_expected < 0 || ref_row_expected >= (PetscInt) n_cells_global) {
    ref_row_expected = 0;
  }

  if (!g_pressure_ctx.ready ||
      g_pressure_ctx.n_cells_local != n_cells_local ||
      g_pressure_ctx.n_cells_global != n_cells_global ||
      g_pressure_ctx.n_faces != n_faces ||
      g_pressure_ctx.ref_row != ref_row_expected) {
    ctx_status = fvm_petsc_pressure_ctx_destroy(&g_pressure_ctx);
    if (ctx_status != 0) {
      return -50 + ctx_status;
    }
    ctx_status = fvm_petsc_pressure_ctx_build(&g_pressure_ctx,
                                              n_cells_local,
                                              n_cells_global,
                                              n_faces,
                                              local_to_global,
                                              face_owner,
                                              face_neighbor_local,
                                              face_neighbor_global,
                                              skip_face,
                                              diagonal,
                                              face_coeff,
                                              reference_cell,
                                              max_iterations,
                                              rtol,
                                              atol);
    if (ctx_status != 0) {
      (void) fvm_petsc_pressure_ctx_destroy(&g_pressure_ctx);
      return ctx_status;
    }
  }

  // Load RHS and initial guess into the PETSc vectors.
   // Fast path (use_direct_solution_copy == 1): compact global IDs are
  if (g_pressure_ctx.use_direct_solution_copy) {
    PetscScalar *rhs_arr = NULL;
    PetscScalar *sol_arr = NULL;

    ierr = VecGetArray(g_pressure_ctx.rhs_vector, &rhs_arr);
    if (ierr != PETSC_SUCCESS) return -15;
    ierr = VecGetArray(g_pressure_ctx.solution_vector, &sol_arr);
    if (ierr != PETSC_SUCCESS) return -25;

    for (row = 0; row < (PetscInt) n_cells_local; ++row) {
      PetscScalar rv = (PetscScalar) rhs[row];
      PetscScalar p0 = (PetscScalar) solution[row];
      if (!isfinite((double) rv)) rv = 0.0;
      if (!isfinite((double) p0)) p0 = 0.0;
      rhs_arr[row] = rv;
      sol_arr[row] = p0;
    }

    ierr = VecRestoreArray(g_pressure_ctx.rhs_vector, &rhs_arr);
    if (ierr != PETSC_SUCCESS) return -15;
    ierr = VecRestoreArray(g_pressure_ctx.solution_vector, &sol_arr);
    if (ierr != PETSC_SUCCESS) return -25;
    /* No VecAssembly needed: values were written directly to owned memory. */
  } else {
    for (row = 0; row < (PetscInt) n_cells_local; ++row) {
      PetscScalar rhs_value = (PetscScalar) rhs[row];
      PetscScalar p0        = (PetscScalar) solution[row];
      if (!isfinite((double) rhs_value)) rhs_value = 0.0;
      if (!isfinite((double) p0))        p0        = 0.0;
      g_pressure_ctx.rhs_values[row] = rhs_value;
      g_pressure_ctx.sol_values[row] = p0;
    }

    ierr = VecSetValues(g_pressure_ctx.rhs_vector, (PetscInt) n_cells_local,
                        g_pressure_ctx.sol_rows, g_pressure_ctx.rhs_values,
                        INSERT_VALUES);
    if (ierr != PETSC_SUCCESS) return -15;
    ierr = VecSetValues(g_pressure_ctx.solution_vector, (PetscInt) n_cells_local,
                        g_pressure_ctx.sol_rows, g_pressure_ctx.sol_values,
                        INSERT_VALUES);
    if (ierr != PETSC_SUCCESS) return -25;
    ierr = VecAssemblyBegin(g_pressure_ctx.rhs_vector);
    if (ierr != PETSC_SUCCESS) return -20;
    ierr = VecAssemblyEnd(g_pressure_ctx.rhs_vector);
    if (ierr != PETSC_SUCCESS) return -21;
    ierr = VecAssemblyBegin(g_pressure_ctx.solution_vector);
    if (ierr != PETSC_SUCCESS) return -25;
    ierr = VecAssemblyEnd(g_pressure_ctx.solution_vector);
    if (ierr != PETSC_SUCCESS) return -25;
  }

  // Detach constant null space when any rank pins a Dirichlet row, else attach; n_dirichlet_rows is fixed per handle so this is cached, not redone every call.
  if (!g_pressure_ctx.global_n_dir_valid) {
    int local_n_dir = (n_dirichlet_rows > 0) ? n_dirichlet_rows : 0;
    int global_n_dir = 0;
    MPI_Allreduce(&local_n_dir, &global_n_dir, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    g_pressure_ctx.global_n_dir = global_n_dir;
    g_pressure_ctx.global_n_dir_valid = 1;
  }
  {
    int want_attached = (g_pressure_ctx.global_n_dir > 0) ? 0 : 1;
    if (g_pressure_ctx.nullsp_attached != want_attached) {
      if (!want_attached) {
        (void) MatSetNullSpace(g_pressure_ctx.matrix, NULL);
      } else {
        MatNullSpace nullsp = NULL;
        PetscErrorCode ierr_ns = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullsp);
        if (ierr_ns == PETSC_SUCCESS) {
          (void) MatSetNullSpace(g_pressure_ctx.matrix, nullsp);
          MatNullSpaceDestroy(&nullsp);
        }
      }
      g_pressure_ctx.nullsp_attached = want_attached;
    }
  }
  (void) dirichlet_rows;
  (void) dirichlet_values;

  /* Null space handles gauge; no ref-row pin needed. */
  ksp_max_iterations = (PetscInt) ((max_iterations > 2000) ? max_iterations : 2000);
  current_rtol = (PetscReal) rtol;
  current_atol = (PetscReal) atol;
  if (!g_pressure_ctx.ksp_tolerances_set ||
      g_pressure_ctx.last_max_iterations != ksp_max_iterations ||
      g_pressure_ctx.last_rtol != current_rtol ||
      g_pressure_ctx.last_atol != current_atol) {
    ierr = KSPSetTolerances(g_pressure_ctx.ksp, current_rtol, current_atol,
                            PETSC_DEFAULT, ksp_max_iterations);
    if (ierr != PETSC_SUCCESS) return -31;
    g_pressure_ctx.last_max_iterations = ksp_max_iterations;
    g_pressure_ctx.last_rtol = current_rtol;
    g_pressure_ctx.last_atol = current_atol;
    g_pressure_ctx.ksp_tolerances_set = 1;
  }

  ierr = KSPSolve(g_pressure_ctx.ksp, g_pressure_ctx.rhs_vector, g_pressure_ctx.solution_vector);
  if (ierr != PETSC_SUCCESS) return -33;
  ierr = KSPGetConvergedReason(g_pressure_ctx.ksp, &reason);
  if (ierr != PETSC_SUCCESS) return -34;
  ierr = KSPGetIterationNumber(g_pressure_ctx.ksp, (PetscInt *) iteration_count);
  if (ierr != PETSC_SUCCESS) return -35;
  {
    Vec true_resid = NULL;
    PetscReal true_resid_norm = 0.0;
    ierr = KSPBuildResidual(g_pressure_ctx.ksp, NULL, NULL, &true_resid);
    if (ierr != PETSC_SUCCESS) return -36;
    ierr = VecNorm(true_resid, NORM_2, &true_resid_norm);
    if (ierr != PETSC_SUCCESS) return -36;
    * residual_norm = (double) true_resid_norm;
    VecDestroy(&true_resid);
  }

  if (reason < 0) {
    PC pc = NULL;
    ierr = KSPGetPC(g_pressure_ctx.ksp, &pc);
    if (ierr != PETSC_SUCCESS) return -29;
    ierr = PCSetType(pc, PCBJACOBI);
    if (ierr != PETSC_SUCCESS) return -30;
    ierr = KSPSetTolerances(g_pressure_ctx.ksp, current_rtol, current_atol,
                PETSC_DEFAULT, 4 * ksp_max_iterations);
    if (ierr != PETSC_SUCCESS) return -31;
    g_pressure_ctx.last_max_iterations = 4 * ksp_max_iterations;
    g_pressure_ctx.last_rtol = current_rtol;
    g_pressure_ctx.last_atol = current_atol;
    g_pressure_ctx.ksp_tolerances_set = 1;
    ierr = KSPSetUp(g_pressure_ctx.ksp);
    if (ierr != PETSC_SUCCESS) return -32;
    /* Zero the initial guess: the failed BoomerAMG solve may have left diverged/NaN
     * values in solution_vector. Starting BJACOBI from x0=0 avoids an immediate
     * KSP_DIVERGED_DTOL at iter=0 caused by r0 = b - A*garbage being enormous. */
    ierr = VecZeroEntries(g_pressure_ctx.solution_vector);
    if (ierr != PETSC_SUCCESS) return -37;
    ierr = KSPSolve(g_pressure_ctx.ksp, g_pressure_ctx.rhs_vector, g_pressure_ctx.solution_vector);
    if (ierr != PETSC_SUCCESS) return -33;
    ierr = KSPGetConvergedReason(g_pressure_ctx.ksp, &reason);
    if (ierr != PETSC_SUCCESS) return -34;
    ierr = KSPGetIterationNumber(g_pressure_ctx.ksp, (PetscInt *) iteration_count);
    if (ierr != PETSC_SUCCESS) return -35;
    {
      Vec true_resid = NULL;
      PetscReal true_resid_norm = 0.0;
      ierr = KSPBuildResidual(g_pressure_ctx.ksp, NULL, NULL, &true_resid);
      if (ierr != PETSC_SUCCESS) return -36;
      ierr = VecNorm(true_resid, NORM_2, &true_resid_norm);
      if (ierr != PETSC_SUCCESS) return -36;
      petsc_residual = true_resid_norm;
      * residual_norm = (double) true_resid_norm;
      VecDestroy(&true_resid);
    }

    if (reason < 0) {
      // Soft failure: accept DIVERGED_ITS with residual within 1e4*rtol or 1e3*atol (near-zero RHS / atol-only noise floor); anything else is fatal.
      if (!(reason == KSP_DIVERGED_ITS && (petsc_residual < (PetscReal)(rtol * 1.0e4) ||
                                            petsc_residual < (PetscReal)(atol * 1.0e3)))) {
        /* Restore the original preconditioner before returning so that the KSP is
         * not permanently left in BJACOBI state, which would cause every subsequent
         * call to fail immediately (all MPI ranks see iter=0 failures in lockstep). */
        if (strncmp(g_pc_solver_type, "boomeramg", 9) == 0) {
          (void) PCSetType(pc, PCHYPRE);
          (void) PCHYPRESetType(pc, "boomeramg");
        } else {
          (void) PCSetType(pc, PCGAMG);
        }
        (void) KSPSetTolerances(g_pressure_ctx.ksp, current_rtol, current_atol,
                                PETSC_DEFAULT, ksp_max_iterations);
        g_pressure_ctx.last_max_iterations = ksp_max_iterations;
        g_pressure_ctx.last_rtol = current_rtol;
        g_pressure_ctx.last_atol = current_atol;
        g_pressure_ctx.ksp_tolerances_set = 1;
        (void) KSPSetUp(g_pressure_ctx.ksp);
        status = -(1000 + (int) reason);
        return status;
      }
    }

    /* Restore the originally configured preconditioner. */
    if (strncmp(g_pc_solver_type, "boomeramg", 9) == 0) {
      ierr = PCSetType(pc, PCHYPRE);
      if (ierr != PETSC_SUCCESS) return -30;
      ierr = PCHYPRESetType(pc, "boomeramg");
      if (ierr != PETSC_SUCCESS) return -30;
    } else {
      ierr = PCSetType(pc, PCGAMG);
      if (ierr != PETSC_SUCCESS) return -30;
    }
    ierr = KSPSetTolerances(g_pressure_ctx.ksp, current_rtol, current_atol,
                            PETSC_DEFAULT, ksp_max_iterations);
    if (ierr != PETSC_SUCCESS) return -31;
    g_pressure_ctx.last_max_iterations = ksp_max_iterations;
    g_pressure_ctx.last_rtol = current_rtol;
    g_pressure_ctx.last_atol = current_atol;
    g_pressure_ctx.ksp_tolerances_set = 1;
    ierr = KSPSetUp(g_pressure_ctx.ksp);
    if (ierr != PETSC_SUCCESS) return -32;
  }

  if (g_pressure_ctx.use_direct_solution_copy) {
    const PetscScalar *sol_arr = NULL;
    ierr = VecGetArrayRead(g_pressure_ctx.solution_vector, &sol_arr);
    if (ierr != PETSC_SUCCESS) return -43;
    for (row = 0; row < (PetscInt) n_cells_local; ++row)
      solution[row] = (double) sol_arr[row];
    ierr = VecRestoreArrayRead(g_pressure_ctx.solution_vector, &sol_arr);
    if (ierr != PETSC_SUCCESS) return -44;
  } else {
    ierr = VecScatterBegin(g_pressure_ctx.sol_scatter, g_pressure_ctx.solution_vector,
                           g_pressure_ctx.sol_local, INSERT_VALUES, SCATTER_FORWARD);
    if (ierr != PETSC_SUCCESS) return -41;
    ierr = VecScatterEnd(g_pressure_ctx.sol_scatter, g_pressure_ctx.solution_vector,
                         g_pressure_ctx.sol_local, INSERT_VALUES, SCATTER_FORWARD);
    if (ierr != PETSC_SUCCESS) return -42;

    {
      const PetscScalar *sol_arr = NULL;
      ierr = VecGetArrayRead(g_pressure_ctx.sol_local, &sol_arr);
      if (ierr != PETSC_SUCCESS) return -43;
      for (row = 0; row < (PetscInt) n_cells_local; ++row) {
        solution[row] = (double) sol_arr[row];
      }
      ierr = VecRestoreArrayRead(g_pressure_ctx.sol_local, &sol_arr);
      if (ierr != PETSC_SUCCESS) return -44;
    }
  }
  return status;
#else
  (void) n_cells_local;
  (void) n_cells_global;
  (void) n_faces;
  (void) local_to_global;
  (void) face_owner;
  (void) face_neighbor_local;
  (void) face_neighbor_global;
  (void) skip_face;
  (void) diagonal;
  (void) rhs;
  (void) face_coeff;
  (void) reference_cell;
  (void) max_iterations;
  (void) rtol;
  (void) atol;
  (void) solution;
  (void) iteration_count;
  (void) residual_norm;
  return -1;
#endif
}
