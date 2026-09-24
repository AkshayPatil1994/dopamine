/* SPDX-FileCopyrightText: 2026 DOPAMINE contributors */
/* SPDX-License-Identifier: AGPL-3.0-only */
/* NVTX range markers for Nsight Systems, called from mod_profiling. NVTX3 is header-only and
 * a no-op unless a profiler injects itself, so leaving these compiled in costs nothing. */
#include <nvtx3/nvToolsExt.h>

void dopamine_nvtx_push(const char *name) { nvtxRangePushA(name); }
void dopamine_nvtx_pop(void) { nvtxRangePop(); }
