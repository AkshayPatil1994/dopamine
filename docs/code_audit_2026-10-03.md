# Code audit — 2026-10-03

Scope: all of `src/` read file-by-file, except `mod_cgns_reader.F90` and the second half of `mod_icem_reader.f90` (offline decompose-tool readers only; both readers have since been removed from `dev`). No code was modified. Three small experiments were run with `build/dopamine`; results were cross-checked against `tests/regression/reference/*.csv`.

Status legend: **[confirmed]** = traced in code (and/or evidence); **[suspected]** = needs a test.

## Status on `dev`

| # | Status | Commit subject |
|---|---|---|
| 1 | fixed (symmetric `d_cf(own)+d_cf(partner)` coefficient on every periodic edge) | Make periodic pressure-matrix coefficients identical on both sides of a seam |
| 2 | fixed as symmetry only (processor faces use the stencil-free reconstruction on both ranks; accuracy at partition seams is first-order-like, a deeper halo would restore it) | Use a stencil-free TVD/QUICK reconstruction on both sides of a processor face |
| 4 | fixed | Apply the n_nonortho_corr auto-upgrades when the input leaves it at auto |
| 5 | fixed (also zero-gradient when the effective diffusivity is 0) | Impose Neumann scalar fluxes through the effective diffusivity |
| 7 | fixed (collective finite check every 10 steps and before output; abort path exercised by injecting a NaN on rank 1: detected at the next check and every rank aborted. Natural blow-ups are caught earlier by the existing time-step-collapse guard) | Abort the run when velocity or pressure turn non-finite |
| 6 | fixed. The stage structure is: predictor `u*` has no pressure, the pre-correction removes `gk_dt*grad(p_old)` using the whole pressure (outlet value `p_val`), and the solve/correction work on the increment phi, whose outlet value must be 0. The code used `p_val` for phi in the solve, the face-flux and velocity corrections, and (for `n_nonortho_corr > 1`) in the boundary values summed into the total pressure (`2 p_val`); a mean-outlet anchor then hid it for equal `p_val`. Now: phi = 0 on outlets for increments, `p_val` only for whole-pressure uses (pre-correction, initial projection), anchor dropped when outlets exist. Hydrostatic outlet box: p = 1 - z exactly; pressure-driven box: Umax = t. The suspected cross-term double count was this bug, not a double count: skewed pressure-driven box Umax was 2x too large for `n_nonortho_corr >= 2`, now 0.125 vs 0.121 | Treat the pressure increment as zero on Dirichlet outlets and drop the mean-outlet anchor |
| 11 | fixed: unknown `patch_type`, scheme/SGS-model/pressure-solver/gradient-scheme names (now case-insensitive), >32 patches, >32 PETSc options and `write_interval = 0` stop with a message; a missing BC patch aborts every rank. Unlisted mesh patches still warn and default to open | Reject unknown patch types...; Match scheme names case-insensitively...; Abort every rank when an input BC patch is missing |
| 14 | fixed (cell-count sum check on read) | Abort when a preprocessed mesh is run on fewer ranks than it was decomposed for |
| 9 | mostly fixed: output/sample files continue and truncate at the restart, missing `.nc` starts a new file, repartitioned restart rebuilds `U_face`, stats CSV trimmed and `E0` kept, the wall-stress filter and the synthetic-eddy state (positions, signs, RNG) are checkpointed, and gradients/boundary values are rebuilt after a restart. Exact (bit-identical) for periodic/wall-model channel restarts; ~5e-6 relative for inflow/outlet/slip cases because the first stage's boundary extrapolation used the previous stage's gradients, which are recomputed rather than stored (storing them would roughly triple the checkpoint). An existing `.nc` with a different cell layout (after repartitioning) is moved to `.nc.old` and a new file is started. Still open: a shorter rerun leaves a stale tail in the unlimited time dimension | Continue output and sample files...; Rebuild the face flux...; Trim the stats CSV...; Checkpoint the wall-stress time filter...; Checkpoint the synthetic-eddy state... |
| 10 | fixed (no skew correction on periodic faces, as in the SoA kernels; latent path, not exercised by the regression cases) | Skip the skewness correction on periodic faces in the host face interpolation |
| 12 | partly fixed: unsupported element types and untagged boundary faces now stop with a message, large arrays moved to the heap. Still open: only the first physical tag per entity, vertex-mean face centres, tet-split volumes | Make the Gmsh reader reject unsupported element types and untagged boundary faces |
| 3 | fixed: the partitioner now ships vertex-adjacent cells across cuts (ghosts with no faces, in the halo lists and as `vx_owned/vx_ghost` pairs, also in `mesh.bin`), and step 1 of every regression case agrees with np=1 to machine precision at np=2,3,4,8 (rb `eps` was 2.5% off). Found on the way: the dynamic-Smagorinsky test filter ignored periodic partners on single-rank meshes (np=1 `nu_sgs_mean` was 13% too low) | Carry vertex-adjacent cells across partition cuts...; Include periodic partners in the dynamic-Smagorinsky test filter... |
| 8 | fixed for the wall-cell imbalance: the pre-correction gradient had no wall-face term (walls, slip, symmetry, inlets), so a wall cell saw half of a hydrostatic gradient; it now adds `A (f.n) n` with `f` = constant/harmonic forcing + Boussinesq buoyancy. Skewed closed box Umax 3.0e-2 -> 5.9e-3; `rb` (Boussinesq, fixed-T walls) no longer gains spurious energy (1.2e-4 -> 5e-5 over 20 steps). What remains in the orthogonal closed box is an O(dt g) startup checkerboard (6.2e-3 at dt = 2e-3, 1.8e-3 at 5e-4) that viscosity does not damp, from the impulsive start with p = 0 | Carry the body force into the wall pressure gradient of the pre-correction |
| 13 | resolved by removing the readers | |

## Other fixes made while working through the list

- `t_end` is now hit exactly (the last step is shortened) and a CFL-clamped `write_time` step is no longer written twice.
- `copy_solver_params` is an intrinsic assignment (the hand-kept field count compared constants and could never fail).
- `mesh_destroy` frees the SoA mirrors and colouring on host builds.
- GPU build audit (OpenACC + CUDA, two RTX A6000): the non-orthogonal cross-flux kernels called a helper taking pointer arrays inside an `acc routine seq` (illegal device address, `pdrop`/`hydro`/`hydrowall`/`inflow*`); the owner weight is now inlined. `Tcell` was missing from the data clause of the wall pre-correction (present-table error at the first projection). With both fixed all regression cases pass at np=1,2 on GPU, `inflow` also with `DOPAMINE_GPU_AWARE_MPI=0`, and tgv/inflow/channel at np=4 on 2 GPUs. Multi-node and timing were not run.
- `noise_seed` makes the initial noise reproducible; radius-0 scalar sources and unreadable source/sampling files now warn.
- The mesh-quality report counts processor faces once; the divergence scratch array is on the heap; the MOST wall heat flux is skipped when T is not solved.
- Left alone on purpose: `ksp_max_iterations = max(p_max_iter, 2000)` (a floor that guards AMG robustness), `foam bc_type_to_id` defaults (only used for patches missing from the input, which already warn), int32 mesh headers, and the O(n^2) sampling set-up.

## Suggested fix order

1. #1 periodic pressure coefficients, 2. #2 ghost-upwind TVD/QUICK, 4. #4 `n_nonortho_corr` auto, 5. #5 Neumann scalar BC, 6. #6 outlet `p_val`, 7. #7 NaN abort; then the rest.

## High severity

### 1. Asymmetric periodic pressure coefficients on partitioned meshes [confirmed]
`src/solver/mod_pressure_solver.f90:287-311` (built at 241-330).
- In partitioned meshes `periodic_neighbor` is always populated, so the `pair_face` branch (line ~253) never runs and `pair_face = 0` for *every* periodic face, same-rank pairs included.
- Coefficient falls into the "cross-rank" branch: `d_per = d_proj_o + periodic_pair_dcf(f_bnd)` = projected own distance + partner's **un-projected** `d_cf`.
- Single-rank runs use `d_proj_o + d_proj_pair` (both projected) → np=1 and np>1 differ.
- Cross-rank seams: each rank computes `area/(dpA + dcfB)` vs `area/(dpB + dcfA)`; unequal unless cells are aligned. Matrix is non-symmetric but flagged `MAT_SPD` (`pressure_backend_petsc.c:370-373`); `correct_face_fluxes` gives `ff_A != -ff_B` → mass not conserved at the seam.
- Invisible on aligned Cartesian meshes. Matches the Coceal np>1 failure in memory.
- Fix idea: carry the partner's projected normal distance (`periodic_pair_dproj`) from the partitioner (global mesh has both faces) and use `d_proj_o + d_proj_partner` in all branches; make coefficient identical on both sides.

### 2. TVD/QUICK reconstruction is rank-asymmetric at processor faces [confirmed]
`src/numerics/mod_advection.f90`: `soa_find_far_upwind_dev` (~2541: `if (up > nc) return` → far_up = -1), `soa_quick_reconstruction_ctx_dev` (~2613: stencil only if `up <= nc`), QUICK grad buffer sliced to `1:nc` (`advection_rhs_vector_recon_soa`). Host versions walk the ghost's partial face list.
- Face seen from rank A has owned upwind cell → full TVD; from rank B the same cell is a ghost → first-order. `phi_f` differs → non-conservative, decomposition-dependent.
- Affects VANLEER_TVD (default scalar/temperature scheme), KOREN, QUICK, for velocity and scalars.
- Fix idea: make the reconstruction a function of global geometry symmetric across the face (e.g. exchange far-upwind cell index + value, or ghost-layer depth 2 for TVD/QUICK stencils, or compute processor-face `phi_f` on one rank and send).

### 3. LS gradients depend on the decomposition [confirmed]
`src/mesh/mod_mesh.f90:586-632` — vertex-neighbour extension builds vertex→cell CSR over owned cells only, so cells at partition cuts lose edge/vertex neighbours across ranks.
- Evidence: `tests/regression/reference/rb_np1.csv` vs `rb_np2.csv`: `eps` differs ~2.5% at step 1 (LS gradients); channel `nu_sgs_mean` differs ~2e-3.
- Fix idea: include ghost cells' vertices (send vertex lists for ghost cells) or drop vertex extension at partition cuts consistently.

## Medium

### 4. `n_nonortho_corr` auto-upgrade is dead code [confirmed]
`src/main.f90:248` tests `n_nonortho_corr == 1`; parser default is 0 (`mod_input_parser.f90:190,264`; `reference.input:215`). `detect_nonortho_skip` (`mod_fv_solver.f90:1900`) then resolves 0→1/2. The 60°/45° and refinement-interface (=4) upgrades never fire for default input. With `diag_interface` the loop cap is 8 but only 2 sweeps run.
- Fix idea: treat 0 as "auto" in main (apply upgrades when `== 0`) and let the solver resolve the remainder afterwards.

### 5. Neumann passive-scalar BC is unit-inconsistent [confirmed]
`src/bc/mod_bc.f90:623-624`: `C_bnd = C(o) + flux*d_cf/max(alpha_mol,SMALL)`. Solver diffuses with `alpha_mol + nu_sgs/Pr_t` (`mod_viscous.f90`). Default `scalar_diffusivity = 0` → `1e30` boundary value; with SGS the imposed flux is off by `alpha_eff/alpha_mol`.
- Fix idea: pass `nu_sgs(o)/Pr_t` into the BC (or impose the flux directly in the viscous boundary term).

### 6. Non-zero outlet `patch_p` applied to the increment [suspected, symptom observed]
Dirichlet `c_b*p_val` is added to the φ-solve (`mod_pressure_solver.f90:403`, `solve` ~727) but `p = p_old + φ` accumulates (`mod_fv_solver.f90:1544-1548`); `correct_face_fluxes` uses `pval - p(o)` with p=φ (2667); `anchor_pressure_to_outlet` shifts to mean outlet p.
- Experiment (periodic x,y box, z0 outlet p=1, zL outlet p=0, `fz=-1`): Umax = g·t (free fall) on orth mesh → pressure drop did nothing.
- Fix idea: use `p_val - p_old_bnd` as the φ Dirichlet value (and matching flux/velocity corrections).
- Follow-up experiment: setting the φ outlet value to 0 (`p_bnd` already equals `p_val`, so the increment is 0) changed nothing: Umax stays g·t. Pressure at t=0.12 is flat ≈0.6 with ≈0.5 in the end cells, i.e. both outlets are pinned to the *mean* of the two `p_val`s by `anchor_pressure_to_outlet` and no drop is built. The outlet face flux is `Ub·Sf` (extrapolated velocity, no pressure term) and the solve never sees the imposed difference, so the fix needs a derivation of the outlet term in the predictor/correction pair, not a value swap. With a single `p_val` (all outlets equal, the `duct` case) the current code is consistent.
- Second experiment: with φ_b = 0 *and* `anchor_pressure_to_outlet` disabled when outlets exist (it re-pins the mean outlet-cell pressure to the mean `p_val`, erasing any difference between outlet patches), pressure stays exactly 0 and Umax is still g·t. So the momentum predictor never sees the imposed `p_val` either; `p_val` only enters through the φ solve. A real fix has to put the outlet boundary pressure into the predictor/cell pressure gradient and the face flux consistently, and then drop the anchor for Dirichlet-outlet runs. Not attempted further; needs a short derivation of the intended incremental scheme at outlets first.

### 7. No NaN/divergence abort [confirmed]
NaN RHS/solution silently zeroed (`pressure_backend_petsc.c:773-775,792-793`; `mod_pressure_solver.f90:708-709`); `solver_compute_dt` falls back to `params%dt` on NaN. A blown-up run writes NaN until `t_end`.
- Fix idea: global `max|U|` finite check each N steps (collective) → `MPI_Abort` with message; make PETSc interface report non-finite input.

### 8. Body force / buoyancy at walls is unbalanced [confirmed experimentally]
Wall pressure BC ignores the body-force normal component (`mod_bc.f90:423-430`, `a_wall` unused). Hydrostatic closed box, `fz=-1`, no-slip, orth mesh: Umax ≈ 6e-3 at t=0.12 (should be ~0). Affects Boussinesq/ABL near ground.

**Follow-up (first attempt, reverted; see the status table for the fix).** Adding `dp/dn = f.n` (constant force plus Boussinesq buoyancy) to the total-pressure wall BC left the orthogonal box at exactly 6.249e-3 (that mesh sets `skip_grad_p`, so `p_bnd` never enters a gradient) and made the skewed box worse (3.0e-2 -> 3.5e-2). The BC is therefore not the cause. Likely cause: `correct_velocity_psolver_dual` builds the cell correction from face differences only and gives walls no contribution, while the wall face flux is fixed at 0, so a wall cell keeps the `dt*f` its predicted velocity gained but its wall face never carries. A fix needs that wall-face term designed together with the predictor (not a BC value); left open.

### 9. Restart gaps [confirmed]
- Sampling files always `nf90_clobber` (`mod_sampling.f90:388`) → restart destroys earlier samples.
- Output writer appends after existing records, no truncation to restart time (`reopen_netcdf_file`); restart without an existing `.nc` aborts.
- `dopamine-checkpoint-partition` writes `header(7)=0` (drops `U_face`) → first RK stage after repartition has no convection.
- SEM eddy state, wall-model time filter, `L_bnd` not checkpointed; stats `E0_ref` re-latched.

### 10. Host `face_interpolate_*` skew correction wrong on periodic faces [confirmed, latent]
`mod_gradient.f90:610,669`: `epsilon` uses raw centroid difference across the seam (domain length). SoA paths override (`phi_f = phi_f_plain` for `f > n_int`); host paths don't. Reachable on CPU builds with CDS when refinement faces exist (`ke_budget_on` forces host path in `advection_rhs_vector`).

### 11. Input robustness [confirmed]
- Unknown `patch_type` → silently `BC_ZEROGRADIENT` (`noslip` ≠ `noSlip` = open boundary), `mod_input_parser.f90:1018`; unlisted mesh patches also default open.
- `write_interval = 0` → `mod(step,0)` SIGFPE (`main.f90:385`).
- `>32` BC patches silently truncated (`MAX_BC_PATCHES`); `n_petsc_options > 32` overruns `petsc_options_tmp`.
- Unknown `sgs_model`, `pressure_solver`, `gradient_scheme` silently default.
- `bc_setup`: input-patch-not-in-mesh `error stop` only on rank 0.

### 12. Gmsh reader [confirmed]
- `gmsh_nodes_per_element` default 4: point elements (type 15) and other unsupported types desync the list-directed read (`mod_gmsh_reader.f90:975`).
- Boundary entity without physical group → `phys_tag = entity_tag` may match wrong group; unmatched gives `face_patch = -1` → `patch_map(-1)` OOB (~545).
- Automatic arrays `vol_indices(n_elem)`, `bnd_indices(n_elem)` on stack (with `-fopenmp` → overflow ~1M elements).
- Only first physical tag per entity used; face centre is vertex mean; hex/prism/pyramid volumes via tet splits (inconsistent with face-based moments on warped cells).

### 13. Parser/main mesh_type mismatch [resolved]
The CGNS and ICEM/Fluent readers (the only reason the parser accepted `icem|fluent|cgns`) were removed from `dev`; the parser now accepts only `gmsh|openfoam|foam|preprocessed`.

### 14. Preprocessed meshes not validated against rank count [confirmed]
`mesh.bin` header lacks nranks; `main.f90:144-155` doesn't check `nranks` vs decomposition.

## Low
- `copy_solver_params` guard (`mod_fv_solver.f90:316-370`) compares constants; cannot fail.
- `mesh_destroy` leaks SoA/coloring arrays and doesn't reset `soa_ready`.
- Loop can overshoot `t_end` (no final-step truncation, float accumulation).
- Fixed-dt CFL clamp breaks `snap_dt_to_write_time` (double write near boundary).
- Scalar `point/line/sphere` sources with default `radius=0` match nothing; missing OBJ file silently ignored.
- Sampling: nearest-centroid owner (no inside-domain test), O(npts·ncells) init with 2 collectives/point, O(n²) dedup/sort of intersections.
- `report_worst_divergence_cell` stack array `div_u(n_cells)`.
- MOST wall heat flux integrated into T even when `solve_temperature=.false.`.
- `ksp_max_iterations = max(p_max_iter, 2000)` ignores smaller user values.
- Noise IC seeded from clock (not reproducible); weather loader prints on every rank.
- `report_mesh_quality` double-counts processor faces in refine counts.
- `partition`: `n_cells_global` etc. stored as int32 in headers (>2^31 cells unsupported).
- `foam bc_type_to_id`: `inletOutlet`→INLET, `empty`→SLIP (only defaults when patch not in input).

## Unverified suspicions (need tests)
- **Double-counted non-orthogonal cross term in the projection.** After pre-correction + Rhie-Chow, `F_pred = F(u*) − dt(ḡprc·k + coeff·Δp_old)` already removes the old cross part; `correct_face_fluxes` then removes `cross(grad p_total)` (jno=1: `grad p_old`) instead of `cross(grad φ)`. Affects `compute_face_fluxes_and_rhs` / `correct_face_fluxes` (`mod_fv_solver.f90:2259-2270, 2415, 2635-2643`). Hydrostatic experiments were confounded by #6/#8. Test: manufactured linear/harmonic-pressure case on a skewed mesh (curvilinear hex generator used: `scratchpad/hs/mk.py`, perturbation `amp*sin(πX)sin(2πY)sin(πZ)` etc., 12³, 32.5° non-orth) checking φ≈0 at steady state.
- Idempotency probe (`debug_diagnostics`, Taylor-Green IC on the 12^3 box, 1 step): relative L2 change 2.0000e-3 on the orthogonal mesh and 2.0014e-3 on the 32.5-degree skewed mesh, for `n_nonortho_corr` 1 and 4. A skew-specific double-counted cross term would show up as a difference between the two; none is visible at this sensitivity, so the suspicion is not supported (nor ruled out). The same probe with a noise IC gives 0.326 vs 0.324, which only shows the collocated projection is not idempotent for non-smooth fields.
- Rotational pressure correction uses full Poisson RHS (incl. RC + cross terms) rather than ∇·u* (`mod_fv_solver.f90:1533`).

## Experiment notes
- Closed no-slip box, `fz=-1`, skew vs orth, `ncorr` 1/4: Umax(step 60) orth 6.2e-3; skew 3.0e-2 (ncorr=1), 2.2e-2 (ncorr=4).
- Outlet-p box (see #6): orth 0.122 for all ncorr; skew 0.154/0.131/0.131 for ncorr 1/2/4.
- Regression refs np1 vs np2 differ by up to 2.5% in `eps` (rb) → see #3.
