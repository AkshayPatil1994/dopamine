---
title: Theory
nav_order: 2
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# DOPAMINE Theory and Implemented Numerics

This document maps the current implementation to the underlying numerical theory,
with emphasis on choices that materially affect robustness, dissipation, and LES
consistency on unstructured meshes.

## Governing Equations and Modeling Scope

DOPAMINE solves the incompressible Boussinesq Navier–Stokes equations with
optional scalar transport, using Einstein summation over repeated indices:

$$
\partial_j U_j = 0
$$

$$
\partial_t U_i + \partial_j (U_j U_i)
= -\frac{1}{\rho}\partial_i p
+ \partial_j \left[ (\nu + \nu_t)(\partial_j U_i + \partial_i U_j) \right]
- \beta (T - T_{ref})\,g_i + F_i^{ext}
$$

$$
\partial_t T + \partial_j (U_j T)
= \partial_j \left[ \left(\alpha_T + \frac{\nu_t}{Pr_t}\right) \partial_j T \right]
$$

$$
\partial_t C + \partial_j (U_j C) = \partial_j \left[ D_C\, \partial_j C \right]
$$

where $$U_i$$ is the velocity vector, $$p$$ is pressure, $$\rho$$ is density, $$\nu$$
is kinematic viscosity, $$\nu_t$$ is the turbulent eddy viscosity, $$\beta$$ is the
thermal expansion coefficient, $$T_{ref}$$ is the reference temperature, $$g_i$$ is
the gravitational acceleration vector, $$F_i^{ext}$$ is any additional external
forcing (constant body force, Coriolis, or harmonic), $$\alpha_T$$ is the thermal
diffusivity, $$Pr_t$$ is the turbulent Prandtl number, $$C$$ is a generic passive
scalar, and $$D_C$$ its (possibly spatially varying) diffusivity.

The spatial discretization is collocated, cell-centered finite volume on general
polyhedral meshes (OpenFOAM/gmsh ingestion path), with distributed-memory halo
exchange for ghost states.

## Temporal Integration and Projection Coupling

### Explicit transport updates

The explicit transport integrator is `TIME_RK3W`: Wray low-storage third-order Runge-Kutta (Williamson 1980) [@williamson1980].

$$
\begin{align*}
  \alpha &= (8/15,\; 5/12,\; 3/4), \\
  \beta  &= (0,\; -17/60,\; -5/12), \\
  \gamma &= \alpha + \beta = (8/15,\; 2/15,\; 1/3).
\end{align*}
$$

Three stages, three pressure projections per step.
Low-memory: only one previous-stage RHS buffer per field; no initial-state backup needed.
Stage time abscissas: $$c = (0,\; 8/15,\; 2/3)$$.

### Incremental fractional-step pressure projection

After each RK3W stage, DOPAMINE performs a pressure
projection in the incremental Chorin/Kim-Moin family [@chorin1968; @kim1985].
The incremental form reads:

$$
\hat{\mathbf{u}} = \mathbf{u}^n + \Delta t\, H(\mathbf{u}^n),
\qquad
\mathbf{u}^* = \hat{\mathbf{u}} - \Delta t\,\nabla p^n,
$$

$$
\sum_{f\in P} a_f\,(\varphi_N - \varphi_P) = -\frac{1}{\Delta t}\sum_{f\in P} q_f^*,
$$

$$
q_f^{n+1} = q_f^* - \Delta t\,a_f\,(\varphi_N - \varphi_P), \qquad
\mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t\,\nabla\varphi, \qquad
p^{n+1} = \varphi.
$$

Subtracting $$\nabla p^n$$ before assembling $$q_f^*$$ makes the Poisson RHS driven
by the divergence of the corrected predictor $$\mathbf{u}^*$$, not the raw
extrapolated velocity $$\hat{\mathbf{u}}$$. The Poisson solve produces the new
pressure $$\varphi$$ directly (not a pressure increment). The no-volume, negative-sign
RHS form matches the PETSc backend assembly convention and keeps the face-flux
correction and Poisson operator algebraically dual [@perot1993].

### Velocity correction and discrete duality

Cell velocities are corrected using the same stencil coefficients used by the
Pressure Poisson operator:

$$
\mathbf{u}_P^{n+1} = \mathbf{u}_P^* - \frac{\Delta t}{V_P}\sum_{f\in P} \varphi_f\,\mathbf{S}_f,
$$

where $$\varphi_f = g_c \varphi_P + (1-g_c) \varphi_N$$ and $$g_c$$ is the same
geometric distance-ratio interpolant used in the Poisson assembly. Discrete
duality between the face-flux correction and cell-velocity correction is
necessary: if these two operators are inconsistent, divergence control can
degrade even when the Poisson residual is small.

### RK3W incremental projection per stage

For each stage $$k = 1, 2, 3$$ of Wray RK3W, the effective sub-step is
$$\delta_k = \gamma_k\,\Delta t$$. The running pressure `self%p` holds
$$p^{k-1/2}$$ from the previous stage's Poisson solve (or $$p^n$$ at $$k=1$$).
The per-stage projection reads:

$$
\mathbf{u}^* = \hat{\mathbf{u}} - \delta_k\,\nabla p^{k-1/2}, \qquad
\nabla^2 \varphi = -\frac{1}{\delta_k}\nabla\cdot\mathbf{u}^*, \qquad
\mathbf{u}^{k+1} = \mathbf{u}^* - \delta_k\,\nabla\varphi.
$$

Because $$\delta_k \ne \Delta t$$ in general, using $$\delta_k$$ (not $$\Delta t$$)
in both the Poisson RHS and the velocity correction is essential for
divergence-free intermediate velocities. After each stage $$p^{k+1/2} \leftarrow \varphi$$
is stored in `self%p` and used as the pre-correction pressure for the next stage.
No explicit blending with $$p^n$$ is required.

## Convective Fluxes and Dissipation Control

### Advection scheme choices

DOPAMINE allows independent advection scheme selection for each transported field:

| Input parameter | Field governed | Struct default | Input-file default |
|-----------------|---------------|---------------|--------------------|
| `advection_scheme` | momentum | `KEP` | `CENTRAL` |
| `temperature_scheme` | temperature $$T$$ | `VANLEER_TVD` | `VANLEER_TVD` |
| `scalar_tvd_scheme` | passive scalars | `VANLEER_TVD` | `VANLEER_TVD` |

For LES with resolved turbulence, `KEP` (kinetic-energy preserving) is
recommended for momentum to minimise numerical dissipation and better preserve
resolved scales [@morinishi1998; @verstappen2003]. `VANLEER_TVD` is the
recommended default for temperature and passive scalars: it is monotone and
preserves scalar boundedness without requiring high-resolution unbounded schemes.

### Implemented advection families

- `CDS`/`CENTRAL`: linear centered interpolation (default; recommended for WMLES).
- `KEP`/`SKEW`: midpoint/skew form for improved discrete kinetic-energy behavior (recommended for LES).
- `QUICK`: geometry-aware quadratic upwind reconstruction.
- `VANLEER_TVD`: MUSCL with van Leer limiter [@vanleer1974; @sweby1984].
- `KOREN_TVD`: MUSCL with Koren limiter [@koren1993; @sweby1984].

### QUICK on unstructured meshes

DOPAMINE QUICK uses a projected streamline coordinate and local least-squares
quadratic fit (reducing to classical QUICK on aligned uniform grids
[@leonard1979]). In the updated implementation, candidate upwind neighbors are
accepted for alignment $$\ge 0.5$$ instead of a much stricter threshold.

Why this matters: on generic polyhedra, very strict alignment rejection makes
QUICK collapse too often toward linear behavior, effectively injecting hidden
diffusion by stencil starvation. The new threshold keeps the method genuinely
third-order-biased in a broader set of practical meshes while retaining
directional relevance.

### TVD far-upwind across periodic faces

For TVD reconstructions, far-upwind lookup now resolves periodic partner cells
when a boundary face has no direct neighbor index. Without this, periodic faces
could silently drop to first-order upwind due to missing stencil completion.
This fix restores limiter behavior expected from Sweby-admissible MUSCL designs
[@sweby1984].

## Gradient Operators and Non-Orthogonal Diffusion

DOPAMINE provides two gradient schemes, selected via `gradient_scheme`:

**Weighted Least-Squares (default, `LEAST_SQUARES`).** For each cell $$c$$,
$$\nabla\phi_c$$ minimises the weighted sum of squared differences to all face
neighbours (owned and ghost) via precomputed stencil weights $$\mathbf{w}_{ck}$$:

$$
\nabla\phi_c = \sum_{k \in \mathcal{S}(c)} \mathbf{w}_{ck}\,(\phi_k - \phi_c).
$$

Stencil weights are assembled once at mesh initialisation and stored in
`msh%ls_weights`. This is the recommended scheme: it is accurate on skewed
meshes, has bounded stencil size, and reuses precomputed structures at zero
gradient-assembly overhead.

**Green-Gauss (`GREEN_GAUSS`).** The gradient is assembled face-by-face:

$$
\nabla\phi_c = \frac{1}{V_c}\sum_{f\in\mathcal{F}(c)} \phi_f\, \mathbf{S}_f,
$$

with an iterative skewness correction when `gg_iterations > 1`. Each additional
sweep adds a gradient-based face correction:

$$
\phi_f = g_c \phi_P + (1-g_c)\phi_N + \langle\overline{\nabla\phi}_f, \boldsymbol{\varepsilon}_f\rangle,
$$

where $$\boldsymbol{\varepsilon}_f = \mathbf{x}_f - (\mathbf{x}_P + g_c\,\mathbf{d}_{PN})$$
is the face skewness vector and $$\overline{\nabla\phi}_f$$ is the linearly
interpolated gradient from the previous sweep. A full halo exchange of
cell-centred gradients is performed after sweep 1 and before each subsequent
sweep so that ghost-cell gradients at MPI rank boundaries are valid for the
skewness correction.

### Non-orthogonal correction in viscous and pressure operators

Both the viscous flux and the pressure Poisson RHS use Jasak's minimum-correction
(deferred cross-diffusion) split [@jasak1996; @moukalled2016]. The face normal
gradient is decomposed as:

$$
(\nabla\phi)_f \cdot \mathbf{n}_f = \underbrace{a_f\,(\phi_N - \phi_P)}_{\text{implicit orthogonal}} + \underbrace{(\nabla\phi)_f^{(n-1)} \cdot \mathbf{k}_f}_{\text{explicit cross-diffusion}},
$$

where $$\mathbf{k}_f = \mathbf{S}_f/|\mathbf{S}_f| - (\mathbf{S}_f\cdot\hat{\mathbf{d}}_{PN})\hat{\mathbf{d}}_{PN}/|\mathbf{S}_f|$$
is the non-orthogonal component of the face area vector. The implicit operator
retains only the orthogonal (two-point) stencil for unconditional diagonal
dominance; the cross-diffusion correction is deferred to the RHS using the
previous iteration/step gradient.

## Pressure Operator on Non-Orthogonal Meshes

### Interior-face coefficient definition

For interior faces, the pressure Laplacian coefficient now uses orthogonal
projected spacing:

$$
a_f = \frac{|\mathbf{S}_f|}{| (\mathbf{x}_N-\mathbf{x}_P)\cdot \mathbf{n}_f |},
$$

with fallback to geometric `d_cf` only when projection degenerates numerically.

Why: using full centroid distance instead of the normal projection weakens
coupling on skewed meshes and introduces a mesh-angle-dependent error in the
pressure operator. The orthogonal projection is the consistent two-point metric
for the normal gradient component [@moukalled2016].

### Periodic pressure coupling

Periodic edges are assembled explicitly, including cross-rank fallback when the
paired face is not locally present. Thickness is then inferred from local
boundary-face distances to avoid full-domain centroid-distance artifacts.
This prevents severe under-coupling of periodic pressure edges in decomposed
runs.

## SGS and Wall-Model Coupling

### Vreman SGS core model

DOPAMINE implements the anisotropic form of the Vreman algebraic eddy viscosity model
[@vreman2004]. The velocity gradient tensor $$\alpha_{ij} = \partial u_j / \partial x_i$$
enters through the filter-tensor-weighted product $$\beta = \alpha^T G \alpha$$, i.e.

$$
\beta_{ij} = \sum_{k=1}^{3}\sum_{l=1}^{3} \alpha_{ki}\, G_{kl}\, \alpha_{lj},
$$

where $$G$$ is the symmetric per-cell filter-width tensor (see below; it reduces to
$$G_{kl} = \delta^2 \delta_{kl}$$, recovering the isotropic Vreman formula, only on a
uniform axis-aligned mesh). The second invariant of $$\beta$$ is:

$$
B_\beta = \beta_{11}\beta_{22} - \beta_{12}^2
          + \beta_{11}\beta_{33} - \beta_{13}^2
          + \beta_{22}\beta_{33} - \beta_{23}^2,
$$

and the SGS viscosity is:

$$
\nu_t = C_v\sqrt{\frac{B_\beta}{\alpha_{ij}\alpha_{ij}}},
$$

with $$C_v = 0.07$$ (`C_VREMAN` in `mod_parameters.f90`). Note that $$G$$ is
absorbed into $$\beta_{ij}$$ — there is no separate scalar $$\Delta^2$$ multiplier
in the formula. $$\nu_t = 0$$ whenever $$\alpha_{ij}\alpha_{ij} = 0$$ or $$B_\beta \leq 0$$.

### Anisotropic filter-width tensor

$$G$$ is built once at solver initialisation per cell $$c$$ (`build_filter_tensor` in
`mod_vreman.f90`) from the area-weighted second moment of each face's offset from the
cell centroid:

$$
G^{raw}_{mn}(c) = \sum_{f \in \mathcal{F}(c)} a_f \, (x_{f,m}-x_{c,m})(x_{f,n}-x_{c,n}),
\qquad G(c) = \left(\frac{V_c^2}{\det G^{raw}(c)}\right)^{1/3} G^{raw}(c),
$$

where $$a_f$$ is the face area and the scale factor normalises $$\det G = V_c^2$$ so
that on an axis-aligned hex $$G$$ reduces to $$\mathrm{diag}(h_1^2, h_2^2, h_3^2)$$ —
the classical per-direction Vreman filter widths. Area weighting (rather than a plain
face count) keeps the tensor invariant to face subdivision (e.g. four coplanar
sub-faces contribute as the one face they replaced), which matters on
hanging-node-refined meshes. A degenerate cell (near-coplanar face offsets, so
$$\det G^{raw}$$ collapses relative to its trace) falls back to the isotropic
$$V_c^{2/3}$$ width. After construction, $$G$$ is smoothed by two passes of an
area-weighted neighbour blend (`smooth_filter_tensor`) so the filter width ramps
smoothly across a refinement interface instead of stepping discontinuously, which
would otherwise show up as a spurious $$\nu_t$$ jump at the interface.

### Dynamic Smagorinsky (`sgs_model = 'DYNAMIC_SMAG'`)

An alternative to Vreman: the Germano-Lilly dynamic procedure
[@germano1991; @lilly1992], implemented in `mod_dynamic_smagorinsky.f90`. The
resolved strain rate $$S_{ij}$$ and its norm $$|S|$$ enter the standard Smagorinsky
form $$\nu_t = C_s^2 \Delta^2 |S|$$, but instead of a fixed $$C_s$$, the Germano
identity is evaluated locally every step:

$$
L_{ij} = \widehat{U_i U_j} - \hat U_i \hat U_j, \qquad
M_{ij} = 2\Delta^2\left(\widehat{|S|S_{ij}} - \alpha^2 |\hat S|\hat S_{ij}\right),
$$

where $$\widehat{(\cdot)}$$ denotes the test filter and $$\alpha = \hat\Delta/\Delta$$
(taken as 2). $$C_s^2$$ follows from Lilly's (1992) least-squares contraction
$$C_s^2 = \langle L_{ij}M_{ij}\rangle / \langle M_{ij}M_{ij}\rangle$$, clipped to
$$\geq 0$$ to rule out backscatter (which a purely dissipative eddy-viscosity closure
cannot represent stably).

The classical formulation averages $$L_{ij}M_{ij}$$ and $$M_{ij}M_{ij}$$ over a
homogeneous direction (e.g. wall-parallel planes in channel flow) before dividing,
which this unstructured solver has no general way to identify. Instead, both the
test filter $$\widehat{(\cdot)}$$ and the $$\langle\cdot\rangle$$ averaging use the
same compact operator: a volume-weighted average of a cell and its immediate
face-neighbours (periodic partners included). This is noisier than true
homogeneous-direction averaging, so treat `DYNAMIC_SMAG` as a locally-clipped,
sanity-tested implementation rather than a fully validated production model —
see `mod_dynamic_smagorinsky.f90`'s module docstring and `tests/test_dynamic_smagorinsky.f90`.

### Gradient reuse by Vreman and viscous terms

Vreman SGS and the viscous stress assembly both consume the cell-centred velocity
gradient $$\alpha_{ij}$$. In `compute_all_gradients` the gradient is computed once
(LS or GG, according to `gradient_scheme`) and stored in `grad_U`. Vreman then
receives this pre-computed gradient via the `grad_U_ext` optional argument,
avoiding a duplicate sweep. This ensures both the SGS model and the viscous
flux react to the same gradient definition near boundaries and across patches.

### Wall-model / SGS dissipation partition

Cells adjacent to `BC_WALL_MODEL` faces are now forced to `nu_sgs = 0`.

Rationale: in equilibrium WMLES, wall stress and SGS stress should not both be
used to model the same unresolved wall-layer dissipation budget. Keeping Vreman
active in first-off-wall model cells can double-count dissipation and damp
resolved near-wall structures. The partition is consistent with standard
wall-modeled LES philosophy [@piomelli2002].

### Wall-velocity sampling strategy

Wall-model input sampling has been simplified to the owner cell only (single
first-off-wall sample), rather than multi-neighbor blending.

Why: blending from neighbors at varying wall-normal heights corrupts the
effective $$(u_\parallel, y_w)$$ pair presented to the log-law inversion,
especially on unstructured anisotropic cells. A single well-defined sample
point is more physically interpretable and reproducible.

## Boundary Conditions and Periodicity Nuances

- Periodic flux reconstruction includes anti-symmetric paired-face reuse where
  possible, with robust fallback reconstruction if only one side is present.
- Wall/shear BC handling enforces dissipative stress orientation against the
  owner-cell tangential velocity to avoid occasional sign inconsistencies in
  long runs.
- Outlet pressure anchoring and PETSc reference-cell logic maintain matrix
  solvability under pure-Neumann-like pressure structures.

## Atmospheric Boundary Layer (ABL) Module

DOPAMINE includes an optional ABL module for realistic simulation of wind-driven
urban and environmental flows with coupled thermal stratification. Three independent
mechanisms enable high-fidelity inflow and forcing:

### Monin-Obukhov Similarity Theory (MOST) Wall Model

The MOST wall model computes surface stress and heat flux based on local
logarithmic velocity and temperature profiles, with stability-dependent
corrections [@monin1954; @businger1971]:

$$
u(z) = \frac{u_*}{\kappa} \left[ \ln\left(\frac{z}{z_{0m}}\right) - \psi_m\left(\frac{z}{L}\right) + \psi_m\left(\frac{z_{0m}}{L}\right) \right]
$$

$$
T(z) = T_s + \frac{T_*}{\kappa} \left[ \ln\left(\frac{z}{z_{0h}}\right) - \psi_h\left(\frac{z}{L}\right) + \psi_h\left(\frac{z_{0h}}{L}\right) \right]
$$

where $$u_*$$ and $$T_*$$ are friction velocity and temperature scale, $$L$$ is the
Obukhov length (stability parameter), and $$\psi_m$$, $$\psi_h$$ are Businger-Dyer
correction functions for momentum and heat that depend on the stratification regime
(stable, neutral, or unstable). The solver dynamically computes $$u_*$$ and $$T_*$$ via
coupled Newton iteration to maintain consistency between the profile and cell values.

The Obukhov length closes the system. Substituting the surface kinematic heat flux
$$\overline{w'T'}_s = -u_* T_*$$ into $$L = -u_*^3 T_{ref} / (\kappa g \, \overline{w'T'}_s)$$
cancels one power of $$u_*$$:

$$
L = \frac{u_*^2 \, T_{ref}}{\kappa \, g \, T_*}
$$

so $$L < 0$$ for unstable stratification (warm surface, $$T_* < 0$$, upward heat flux),
$$L > 0$$ for stable, and $$|L| \to \infty$$ as $$T_* \to 0$$. The neutral limit is
represented by a finite clip rather than infinity, which leaves a small residual
$$\psi$$ contribution. Because the Newton iteration differentiates through
$$L(u_*, T_*)$$, the momentum and heat balances are solved as a genuinely coupled
system rather than by lagging the stability correction.

$$L$$ is stored per boundary face, not as a domain scalar. Over heterogeneous
terrain the surface fluxes differ enough between adjacent land covers that a
single averaged $$L$$ is not meaningful — a lake would otherwise impose its own
stability on neighbouring vegetation. Each face reuses its own converged $$L$$
as the next step's initial guess; the domain mean is retained only as a
diagnostic. Roughness lengths $$z_{0m}$$, $$z_{0h}$$ are likewise per patch, so
water, vegetation, and built terrain can carry their own values within one run.
Where $$z_{0h}$$ is not prescribed it follows from $$z_{0m}$$ via the inverse
Stanton number, $$\mathrm{k}B^{-1} = \ln(z_{0m}/z_{0h}) \approx 2$$.

**Advantages:**
- Single-layer model: no need for synthetic wall-resolved cells
- Full thermodynamic coupling: temperature-dependent stability corrections
- Dynamic Obukhov length: self-adjusts to local stratification state, per face
- Heterogeneous surfaces: per-patch roughness for mixed water/vegetation/terrain
- Includes all stability regimes with empirically validated Businger-Dyer functions

### Synthetic Eddy Method (SEM) for Turbulent Inflow

The inflow generator is a divergence-free synthetic eddy method in the spirit of
[@poletto2013]. Each eddy is the *curl* of a compact vector potential, so its
velocity is solenoidal by construction rather than by correction:

$$
\mathbf{u}'(\mathbf{x}) = \sum_{e=1}^{N_e} \nabla\phi\!\left(s_e\right) \times \mathbf{a}_e,
\qquad s_e^2 = \sum_{i}\left(\frac{x_i - x_{e,i}}{\sigma_{e,i}}\right)^2,
\qquad \phi(s) = (1-s^2)^3 \ (s<1)
$$

$$\nabla\cdot(\nabla\phi\times\mathbf{a}) = 0$$ identically, for *any* semi-axes
$$\sigma_{e,i}$$ and any strength $$\mathbf{a}_e$$. That is what allows every eddy to carry
its own height-dependent size and stress-fitted strength without losing
$$\nabla\cdot\mathbf{u}'=0$$. The kernel is $$C^1$$, so vorticity is continuous. (Rescaling a
solenoidal field by a Cholesky factor, as in the classical SEM, destroys this property.)

**Stresses.** The strengths are $$\mathbf{a}_e = C\boldsymbol{\xi}$$ with $$\xi_i=\pm1$$ random and
$$CC^T = S$$. Because the kernel is even, the Reynolds stresses delivered by a population of
eddies with density $$N_\sigma/(\sigma_x\sigma_y\sigma_z)$$ follow in closed form:

$$
\overline{u'u'} = Q_y S_{zz} + Q_z S_{yy},\quad
\overline{v'v'} = Q_z S_{xx} + Q_x S_{zz},\quad
\overline{w'w'} = Q_y S_{xx} + Q_x S_{yy},
$$

$$
\overline{u'v'} = -Q_z S_{xy},\quad \overline{u'w'} = -Q_y S_{xz},\quad \overline{v'w'} = -Q_x S_{yz},
\qquad Q_i = \frac{N_\sigma Q_J}{\sigma_i^2}
$$

with $$Q_J = 48\pi\cdot\tfrac{384}{45045}$$. At every eddy's height the solver inverts this
$$3\times3$$ system for $$S$$ from the local target stresses and kernel scales. Not every
(stress, scale) pair is attainable with divergence-free eddies: $$S$$ must be positive
semi-definite, which for isotropic eddies requires $$\overline{u'u'}\le\overline{v'v'}+\overline{w'w'}$$.
Neutral-ABL anisotropy ($$\sigma_u:\sigma_v:\sigma_w=2.4:1.9:1.25$$) is attainable once the eddies
are elongated, which the default scale ratios $$L_x:L_y:L_z = 1:0.5:0.35$$ do. If a target is not
attainable the nearest realisable $$S$$ is used and the setup prints the worst stress error.

**Length scales.** Each eddy's semi-axes are $$\sigma_i = 1.678\,L_i(z)$$, where $$L_i$$ are
integral length scales (the factor makes $$L_x$$ the streamwise integral scale of
$$u'$$ for the isotropic kernel; for strongly elongated eddies the realised
scale differs by tens of percent, see the unit test). The scales can follow a table or
the power law $$L_x(z)=L_{ref}(z/z_{ref})^m$$, with $$L_y,L_z$$ fixed fractions of $$L_x$$. Eddy
density is $$\propto 1/(\sigma_x\sigma_y\sigma_z)$$ so every height holds the same number of eddies per
kernel, which keeps the stresses unbiased where the scales change.

**Convection.** Eddies move along the inlet normal at their own mean speed $$U(z_e)$$
(Taylor's hypothesis with shear). Eddies leaving the recycle box are re-inserted with heights
drawn flux-weighted by $$U(z)$$, which keeps the height distribution stationary.

**Inputs.** Targets come from, in order of precedence: a table (`sem_profile_file`), the
wind-tunnel style laws in the namelist (mean profile + turbulence intensities), or the
Panofsky–Dutton neutral surface layer ($$\sigma_{u,v,w}/u_*=2.4,1.9,1.25$$,
$$\overline{u'w'}=-u_*^2$$).

#### Pre-computed inflow planes

For a pre-computed database the field is sampled on a spanwise-periodic $$(y,z)$$ grid at frame
spacing $$\Delta t_p$$ and replayed with Taylor's hypothesis, i.e. the plane is the inlet. Replay
period $$m$$ uses a randomly chosen stored segment $$k_m$$ shifted by $$s_m$$ nodes,

$$
\mathbf{u}'(y,z,t) = \sin\theta\;\mathbf{u}_{k_m}(y+s_m\Delta y,\,z,\,t-mP) + \cos\theta\;\mathbf{u}_{k_{m-1}}(y+s_{m-1}\Delta y,\,z,\,t-mP+P),
\quad \theta=\frac{\pi}{2}\frac{t-mP}{T_x}\ (t-mP<T_x),
$$

with $$P = T_{seg}-T_x$$ and $$s_m=\lfloor \mathrm{frac}(\varphi_0 + m/\phi)\,N_y\rfloor$$ ($$\phi$$ the golden ratio).
The shifts spread structures over the span; they cannot move the spanwise-mean part. If segment $$k$$
has spanwise-mean waveform $$\bar u_k(z,\tau)$$, the replay contains the periodic part
$$\frac1K\sum_k \bar u_k(z,\tau)$$, whose power is about $$1/K$$ of that of a single segment and is
concentrated at $$1/P$$. Removing it by $$\bar u_k \to \sqrt{K/(K-1)}\,(\bar u_k - \frac1K\sum_j\bar u_j)$$
leaves the variance unchanged and the periodic part exactly zero.

Interpolation from the stored lattice to the inlet faces is a linear low-pass operator. Variance
retained on a wave with $$n$$ points per wavelength, averaged over offsets:

| $$n$$ | linear | cubic | 6-point Lagrange |
|---|---|---|---|
| 4 | 67% | 89% | 94% |
| 6 | 83% | 97% | 99% |
| 8 | 90% | 99.1% | 99.9% |

so linear interpolation is not used. In the spanwise direction the nodes are periodic and
uniform, where trigonometric interpolation is exact for the resolved wavenumbers; faces that
coincide with stored nodes are copied unchanged.

### Davies Relaxation (Nudging)

Large-scale atmospheric forcing is applied via Newtonian relaxation
to time-varying reference profiles [@davies1976]:

$$
\left.\frac{\partial \mathbf{U}}{\partial t}\right|_{nudge} = -\frac{1}{\tau}(\mathbf{U}_{cell} - \mathbf{U}_{ref}(z,t))
$$

$$
\left.\frac{\partial T}{\partial t}\right|_{nudge} = -\frac{1}{\tau}(T_{cell} - T_{ref}(z,t))
$$

where $$\tau$$ is the relaxation timescale (typically 1800 s for weather-model
coupling) and $$\mathbf{U}_{ref}$$, $$T_{ref}$$ are interpolated from a loaded
profile file (CSV or netCDF format). Nudging is applied only above a
user-specified height $$z_{start}$$, permitting localized LES evolution
in the lower layers while coupling to synoptic-scale motion above.

**Integration:**
- MOST fluxes are computed once per RK stage and used directly as source terms
- SEM fluctuations are added to inlet boundary conditions
- Nudging source is integrated into the explicit RHS at each step
- All features are optional and gated by input flags; backward compatibility
  is preserved for non-ABL simulations

## Mesh Partitioning for Parallel Runs

For MPI-parallel runs, DOPAMINE uses a two-level partitioning strategy to
balance computational load while minimizing memory and communication overhead.

### Default: space-filling curve partitioning

The partitioner orders cells along a 3D space-filling curve and cuts that
ordering into equal pieces [@bader2012]. An SFC gives a locality-preserving
ordering, so geometrically nearby cells land on nearby ranks and the number of
cross-rank face neighbours stays small.

The curve used is Morton (Z-order): at each level the octant index is
$$(r_x \ll 2)\,|\,(r_y \ll 1)\,|\,r_z$$ with no per-level rotation state,
which makes the key cheap to compute and radix-sortable. Measured against a true
Hilbert construction on periodic structured grids this makes almost no
difference: the two produce identical cuts on power-of-two cubes and land within
1–3% of each other otherwise, with Morton often yielding slightly fewer
neighbours per rank. Both sit within about 6% of the theoretical minimum surface
area of ideal Cartesian blocks, so partition shape is not a limiting factor for
box-like domains.

Cells are assigned so that each rank owns an equal *count*, since solver cost
scales with cell count rather than volume. Equalising volume instead is
available as an option but imbalances graded meshes badly: on a tanh
wall-graded mesh the near-wall ranks end up with roughly 3.5 times the mean
cell count at 1024 ranks, and a bulk-synchronous step is paced by its slowest
rank.

**Advantages:**
- Single-pass, non-iterative algorithm with O(n log n) time complexity
- Memory-efficient: no global graph construction, unlike multilevel graph partitioners
- Scalable to very large rank counts (tested up to >1,000 ranks)
- Load-balanced by cell count, with near-minimal surface area on box domains

## PETSc Backend and Linear Solve Layer

Pressure solves are delegated to the PETSc-backed sparse operator path
[@balay2021petsc]. The Fortran side constructs local-to-global maps, diagonal
and face-coupling coefficients, and skip masks for invalid/degenerate faces;
the backend then performs Krylov/PC solve using these assembled arrays.

## Verification: Spatial Convergence

The convergence suite in `examples/convergence_test/` measures the observed order of
accuracy on the 3D Taylor-Green vortex (periodic $$[0,2\pi]^3$$, $$\nu = 0.01$$, plain
`CENTRAL` advection, `RK3W`, no SGS model, pressure tolerance $$10^{-10}$$). The raw numbers
are in [results.json]({{ site.baseurl }}/assets/convergence/results.json). To reproduce:

```bash
cd examples/convergence_test
python3 run_convergence.py all --bin ../../build/dopamine
python3 plot_convergence.py results/results.json --out ../../docs/assets/convergence
```

There is no closed-form solution of the 3D problem for $$t > 0$$, so two kinds of
measurement are used: self-convergence on successive meshes, which needs no
exact solution, and the analytic initial values $$E(0) = 1/8$$ and
$$\varepsilon(0) = 0.75\,\nu$$. Temporal accuracy is discussed separately below.

### Spatial order: second order

![Spatial convergence]({{ site.baseurl }}/assets/convergence/spatial_convergence.png)

Meshes are $$N^3$$ hexahedra with $$N = 8, 16, 32, 64$$, one fixed $$\Delta t = 4\times10^{-3}$$
for all of them, so the time error is the same on every mesh and drops out of the
differences. The `warped` family is the same mesh pushed through a smooth periodic
coordinate map (maximum non-orthogonality about 10 degrees), which exercises the
non-orthogonal corrections and least-squares gradients; its exact solution is unchanged.

| Quantity | Family | Observed order (successive refinements) |
|---|---|---|
| $$\lvert E_N - E_{2N}\rvert$$ at $$t=1$$ | hex | 1.74, 2.13 |
| $$\lvert E_N - E_{2N}\rvert$$ at $$t=1$$ | warped | 1.86, 2.17 |
| $$\lvert\varepsilon_0 - 0.75\nu\rvert$$ (true error) | hex | 1.52, 1.87, 1.97 |
| $$\lvert\varepsilon_0 - 0.75\nu\rvert$$ (true error) | warped | 1.52, 1.87, 1.97 |

Both families approach the design second order as $$h \to 0$$. The coarsest meshes are
pre-asymptotic (order 1.5 at $$N = 8 \to 16$$).

### Discussion: temporal accuracy

The Runge-Kutta stages are third order (shown with the model below), but the solver does not show third order in $$\Delta t$$ at
a fixed mesh: the observed order is first, with an error that scales as $$\Delta t\,h^2$$. The
measurements and their cause are given here. The corresponding
[figure]({{ site.baseurl }}/assets/convergence/temporal_convergence.png) is produced by
`plot_convergence.py`.

Time steps are halved from $$0.2$$ to $$0.0016$$ at fixed $$t_{end} = 1$$ on $$16^3$$, $$32^3$$ and
$$64^3$$ meshes (largest step limited to $$U\Delta t/h \le 0.9$$).

| Mesh | Observed order of $$\lvert E(\Delta t) - E(\Delta t/2)\rvert$$ |
|---|---|
| $$16^3$$ | 0.89, 0.95, 0.97, 0.99, 0.99, 1.00 |
| $$32^3$$ | 0.95, 0.97, 0.99, 0.99, 1.00 |
| $$64^3$$ | 0.98, 0.99, 1.00, 1.00 |

**The observed temporal order is 1, not the 3 of the Wray scheme.** This holds from
$$\Delta t = 0.2$$ down to $$0.003$$ and is insensitive to the advection scheme, the
rotational-projection switch, the gradient scheme and the pressure tolerance (checked on the
$$16^3$$ mesh). In an inviscid run ($$\nu = 10^{-8}$$, panel b), where the exact answer is
$$E = 1/8$$ for all time, the solver loses energy at a rate $$c(h)\,\Delta t$$, and $$c(h)$$
falls by a factor of 4.4 and 4.1 for each halving of $$h$$ (panel c): the dominant time-step
error is $$O(\Delta t\,h^2)$$. The scheme is consistent (the term vanishes under joint
refinement), but refining $$\Delta t$$ alone on a fixed mesh gains only a factor of 2 per halving.

#### Source of the $$O(\Delta t\,h^2)$$ term

`examples/convergence_test/projection_model.py` is a small numpy model of one RK3W step
(cell velocities plus a separate face flux, the stage weights above, wide-gradient
pre-correction, Rhie-Chow face flux, compact Poisson solve, wide-gradient velocity
correction). It reproduces the solver's inviscid energy loss to 0.03% (5.024e-6 against
5.025e-6 at $$N=16$$, $$\Delta t=0.01$$), and the same ablation run in the solver gives the same
numbers. The findings:

1. The face flux is $$q^* = \mathrm{interp}(\hat{\mathbf u}) - \delta\,\nabla_c p^{old}$$ (compact
   gradient) because the Rhie-Chow term adds back exactly the cell gradient that the
   pre-correction removed. The Poisson solve therefore does not depend on $$p^{old}$$ at all
   (switching the pre-correction off changes nothing), and the projection of the flux is
   exact for the compact operator $$L_c$$.
2. The cell velocity is corrected with the *wide* (cell-centred) gradient $$G_w$$ of the same
   $$\varphi$$. The adjoint of $$G_w$$ is minus the interpolate-then-difference operator, not
   $$L_c$$, so the cell-velocity update is not an orthogonal projection. Per stage it changes the
   energy by $$-\delta_k^2(a - b/2)$$ with $$a = \langle|\nabla_c p|^2\rangle$$ and
   $$b = \langle|G_w p|^2\rangle$$, where the consistent value is $$-\delta_k^2 b/2$$. The
   difference $$a - b = O(h^2)$$ is dissipated, giving a loss rate of about
   $$\tfrac12\sum_k\gamma_k^2\,\Delta t\,(a-b)$$ (within 10% of the measured value on the
   $$16^3$$ and $$32^3$$ meshes).

| Variant (inviscid TGV, $$t=0.2$$, $$N=16$$) | Energy loss, $$\Delta t = 0.01$$ | Order in $$\Delta t$$ |
|---|---|---|
| as implemented | 5.0e-6 | 1.0 |
| Rhie-Chow term removed | 3.6e-8 (also 3.7e-8 in the solver) | 2.0 |
| cell velocity projected orthogonally | 4.0e-11 | 3.1 |

Making the cell-velocity projection orthogonal restores third order in $$\Delta t$$ in the
model, so the Runge-Kutta stages themselves are third order and the limit is the
projection. Removing the Rhie-Chow term has the same effect but is not an option, since it
suppresses pressure-velocity decoupling. The model is a uniform Cartesian grid; on general
meshes the same mismatch exists but its size was not measured.

## Verification: Advection Schemes

The scalar advection schemes (`temperature_scheme`: `CDS`, `QUICK`, `KOREN_TVD`, `VANLEER_TVD`)
are compared on a periodic bar of unit length with uniform velocity $$U = 1$$, so the exact
solution after one revolution ($$t = 1$$) is the initial profile. Diffusion and SGS are off, the
temperature is a passive scalar (`boussinesq = .true.`, `beta = 0`), time stepping is `RK3W` at
CFL 0.2, and the mesh is $$N \times 4 \times 4$$ cubic hexahedra. The initial profile is a Gaussian
(centre $$0.25$$, $$\sigma = 0.05$$), optionally with a square pulse on $$[0.55, 0.75]$$
(`ic_type = 'scalar_gauss'` and `'scalar_pulse'`). To reproduce:

```bash
cd examples/convergence_test
python3 run_advection.py --bin ../../build/dopamine
python3 plot_advection.py results/advection.json --out ../../docs/assets/convergence
```

![Advection schemes]({{ site.baseurl }}/assets/convergence/advection_schemes.png)

With 100 cells the schemes behave as their construction predicts. The raw numbers are in
[advection.json]({{ site.baseurl }}/assets/convergence/advection.json).

| Scheme | Final min / max (exact: 0 / 1 for the square pulse, Gaussian peak 1) | Character |
|---|---|---|
| `CDS` | -0.27 / 1.32 | no dissipation; dispersive oscillations around both features |
| `QUICK` | -0.14 / 1.14 | much less dispersion, but over- and undershoots at the pulse edges |
| `KOREN_TVD` | 0.00 / 1.00 | bounded, no oscillations; smears the pulse and clips the Gaussian peak to 0.87 |
| `VANLEER_TVD` | 0.00 / 1.00 | bounded; slightly more smearing than Koren |

![Advection convergence]({{ site.baseurl }}/assets/convergence/advection_convergence.png)

For the smooth Gaussian ($$N = 25$$ to $$400$$) the observed orders over the last refinement are:

| Scheme | $$L_2$$ order | $$L_\infty$$ order |
|---|---|---|
| `CDS` | 1.99 | 2.03 |
| `QUICK` | 2.09 | 2.14 |
| `KOREN_TVD` | 1.61 | 1.32 |
| `VANLEER_TVD` | 1.73 | 1.42 |

`CDS` is second order as designed. `QUICK` also converges at second order here, not the third
order of its one-dimensional finite-difference form. The TVD schemes are between first and
second order, which is expected: the limiter reverts to first order at the Gaussian's extremum,
and this loss is largest in the maximum norm. These results are for the scalar path only; the
momentum advection schemes use separate code and are covered by the Taylor-Green study above.

## Practical Interpretation for Users

For current defaults and updated implementation:

1. Input-file default advection is `CENTRAL` (most robust); struct default is `KEP`.
   Use `KEP` explicitly for resolved LES. Temperature and passive scalars default
   to `VANLEER_TVD` and are controlled independently via `temperature_scheme` and
   `scalar_tvd_scheme`.
2. Dissipation responsibility is explicitly partitioned between SGS and
   wall-model closures near walls.
3. Non-orthogonality handling uses Jasak minimum-correction / deferred cross-diffusion
   in both viscous momentum and pressure operators [@jasak1996].
4. Periodic coupling uses geometric distance-ratio $$g_c$$ rather than a hardcoded
   $$0.5$$ interpolant on non-uniform meshes.
5. Parallel mesh partitioning uses a Morton (Z-order) SFC by default for scalability.

These choices are aimed at robust long-time LES behavior on realistic
unstructured meshes where hidden discretization inconsistencies can dominate
accuracy more than formal order alone.
