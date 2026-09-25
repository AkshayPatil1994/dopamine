/* SPDX-FileCopyrightText: 2026 DOPAMINE contributors */
/* SPDX-License-Identifier: AGPL-3.0-only */
/* Whether MPI can be handed device pointers (CUDA-aware MPI), queried by mod_device.
 * DOPAMINE_GPU_AWARE_MPI=0/1 overrides; otherwise Open MPI's MPIX_Query_cuda_support() decides
 * (compile-time support is not enough: HPC-X reports it at runtime, e.g. when UCX lacks CUDA).
 * Any other MPI defaults to 0, i.e. halos are staged through the host. */
#include <mpi.h>
#include <stdlib.h>
#if defined(OPEN_MPI) && OPEN_MPI
#include <mpi-ext.h>
#endif

int dopamine_mpi_gpu_aware(void) {
  const char *env = getenv("DOPAMINE_GPU_AWARE_MPI");
  if (env != NULL && env[0] != '\0') return env[0] != '0';
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  return MPIX_Query_cuda_support();
#else
  return 0;
#endif
}
