# Test-backed Examples

```@meta
CurrentModule = OWENSFEA
```

OWENSFEA.jl currently keeps its maintained examples in `test/` rather than in a
separate `examples/` directory. The examples below are extracted from those
tests and are the safest starting points for new documentation or regression
cases.

## Straight Cantilever Modal Analysis

Source: `test/CantileverBeamModal.jl`

This case builds a 0.5 m rectangular steel cantilever with 40 linear beam
elements. The root node is fixed in all six DOFs, gravity is disabled, and the
model runs [`modal`](@ref) with `analysisType = "M"`.

The test checks three levels of behavior:

- OWENSFEA frequencies are compared against GXBeam frequencies within 1 percent
  for the first five paired modes.
- The first three two-dimensional bending modes are compared against the
  Euler-Bernoulli formula `sqrt(EI/(rho*A)) / (2*pi) * (k/L)^2`, where
  `k = [1.875, 4.694, 7.855]`.
- Setting `returnDynMatrices = true` writes a MAT file containing assembled and
  boundary-condition-reduced stiffness, mass, and damping matrices.

Minimal modal call pattern:

```julia
feamodel = OWENSFEA.FEAModel(;
    analysisType = "M",
    dataOutputFilename = "none",
    joint,
    gravityOn = false,
    pBC,
    numNodes = mesh.numNodes,
)

freq, damp = OWENSFEA.modal(feamodel, mesh, el)[1:2]
```

When comparing mode lists, remember that the state-space eigenvalues are paired.
Several tests compare every other OWENSFEA frequency.

## Static Cantilever Deflection And Strain

Source: `test/CantileverBeamDisplacement.jl`

This case reuses the 0.5 m cantilever and sweeps tip loads from `1e4` N to
`1e5` N. It runs both linear and nonlinear OWENSFEA paths and compares against
GXBeam plus analytical beam formulas.

The linear analytical checks are:

```julia
tip_displacement = P * L^3 / (3 * E * Iyy)
top_strain = P * x_from_tip * h / (2 * E * Iyy)
```

The test stores all strain channels from [`ElStrain`](@ref):

- `epsilon_x`, `epsilon_y`, and `epsilon_z`;
- `kappa_x`, `kappa_y`, and `kappa_z`.

It compares `kappa_y * h / 2` against the top-fiber bending strain. This is a
useful regression pattern when changing sign conventions or local element
frames.

## Transient Tip-load Response

Source: `test/UnsteadyBeam.jl`

The transient validation case uses a 60 m vertical cantilever, Newmark-beta
analysis (`analysisType = "TNB"`), Rayleigh damping coefficients of `0.005`, and
a sinusoidal tip force. Each step creates [`DispData`](@ref), runs
[`structuralDynamicsTransient`](@ref), and feeds [`DispOut`](@ref) state back
into the next step.

The test compares root force and moment channels against GXBeam time histories.
The current tolerances are intentionally broad and should be treated as
regression diagnostics, not close physical validation.

Skeleton step:

```julia
disp_data = OWENSFEA.DispData(u_s, udot_s, uddot_s, u_sm1)

el_strain, disp_out, reaction = OWENSFEA.structuralDynamicsTransient(
    feamodel,
    mesh,
    el,
    disp_data,
    Omega,
    OmegaDot,
    time,
    delta_t,
    el_storage,
    Fexternal,
    Fdof,
    CN2H,
    rbData,
)

u_s = disp_out.displ_sp1
udot_s = disp_out.displdot_sp1
uddot_s = disp_out.displddot_sp1
```

## Swept Rotating Beam

Sources: `test/CantileverBeamRotating.jl` and
`test/CantileverBeamRotatingModal.jl`

These tests build a two-segment beam with a 45 degree sweep in the second
segment. They exercise orientation, joint welding between beam segments,
rotating-frame loads, nonlinear static spin-up, and Campbell-diagram mode
tracking.

Important details:

- `mesh_beam` in `test/testdeps.jl` creates the welded two-segment mesh and the
  joint row.
- Rotating modal tests set `spinUpOn = true` and `nlOn = true` before calling
  [`autoCampbellDiagram`](@ref).
- `rotSpdArrayRPM = collect(0:1000:6000)` is converted to Hz for OWENSFEA and to
  rad/s for the GXBeam branch.
- The rotating modal test passes explicit `GAy` and `GAz` values so OWENSFEA and
  GXBeam share the same transverse shear stiffness.

Use this case when changing orientation transforms, spin terms, or
`rotationalEffects`.

## Reduced-order Transient Dynamics

Source: `test/NonlinearROM.jl`

The ROM test builds a straight aluminum cantilever, runs full nonlinear
Newmark-beta dynamics, builds a reduced model with [`reducedOrderModel`](@ref),
and compares [`structuralDynamicsTransientROM`](@ref) against the full nonlinear
transient result.

Pinned behavior:

- the nonlinear full-order response differs measurably from the linear
  full-order response;
- the nonlinear ROM stays within 5 percent relative norm of the nonlinear
  full-order tip history;
- the final ROM tip displacement is within 5 percent relative tolerance and
  `5e-5` absolute tolerance.

This is the best maintained example for updating `DispData` with modal states:
when using ROM, include `eta_s`, `etadot_s`, and `etaddot_s` in the next
[`DispData`](@ref).

## Utility Contracts

Source: `test/HelperFunctions.jl`

The helper tests are small enough to read before editing utility code. They pin:

- Gauss quadrature points and weights from [`getGP`](@ref);
- linear and quadratic shape functions from [`calculateShapeFunctions`](@ref);
- type normalization in `Mesh` and [`Ort`](@ref);
- default and explicit transverse shear stiffness;
- assembly dimension checks;
- joint dependent/active DOF maps;
- concentrated nodal term parsing;
- static boundary-condition application;
- modal frequency extraction and reduced displacement reconstruction;
- nonlinear selective-stiffness guardrails.

If a future example needs to explain a low-level map or matrix operation, copy
the smallest relevant pattern from this file rather than relying on a large
coupled turbine case.
