# Theory, Frames, and Units

```@meta
CurrentModule = OWENSFEA
```

## Element Mechanics

OWENSFEA uses a Timoshenko beam formulation. Each node has six degrees of
freedom:

```text
ux, uy, uz, theta_x, theta_y, theta_z
```

The element computes axial, bending, torsional, and transverse shear terms, plus
mass, damping, gravity/body-load, spin, and concentrated-nodal contributions
when those options are active. Element stiffness, mass, damping, and force terms
are evaluated in local beam coordinates and assembled into the global structural
system.

The strain recovery path returns six channels per element through
[`ElStrain`](@ref): `epsilon_x`, `epsilon_y`, `epsilon_z`, `kappa_x`,
`kappa_y`, and `kappa_z`. The cantilever tests compare bending strain with
`kappa_y * h / 2`, so changes to curvature signs should update tests and docs
together.

## Frames

The mesh coordinates define the global structural frame used by public
displacement and load vectors. Element orientation is carried by [`El`](@ref)
fields:

- `psi`: sweep angle about local/global 3-style orientation logic;
- `theta`: cone/elevation angle;
- `roll`: element roll or twist used by the element transform.

The test helper `calculateElementOrientation` computes `Psi_d` and `Theta_d`
from the vector from node 1 to node 2 and stores orientation angles in degrees.
[`SectionPropsArray`](@ref) `twist` values are passed separately and are treated
as section-property data, not mesh orientation metadata.

Local section offsets such as `ycm` and `zcm` are local-section center-of-mass
offsets. The spinning-beam diagnostics verify that steady rotating-frame
inertial loads are evaluated at the center of mass so offset forces and moments
balance correctly about the spin axis for the pinned cases.

`CN2H` is the inertial-to-hub transform used by static and transient calls for
gravity/body loads. The vector `rbData` uses:

```text
1:3  translational acceleration
4:6  angular velocity
7:9  angular acceleration
```

The angular entries in `rbData` are in rad/s and rad/s^2.

## Units

Use SI units unless a file format explicitly states otherwise:

| Quantity | Unit |
| --- | --- |
| length, displacement | m |
| velocity | m/s |
| acceleration | m/s^2 |
| force | N |
| moment | N m |
| mass per unit length | kg/m |
| internal angular speed | rad/s |
| scalar `Omega` and `OmegaDot` inputs | rev/s and rev/s^2 |
| bending stiffness | N m^2 |
| torsional stiffness | N m^2 |
| axial stiffness | N |
| orientation angles in `Ort`/`El` | degrees |
| section-property twist | radians |

## Rotations

Internal displacement rotations and rigid-body angular velocities are in
radians. Mesh orientation helpers and `Ort` fields use degrees for readability
and legacy input compatibility. Convert at the IO boundary and keep internal
solver state in radians.

The public scalar `Omega` and `OmegaDot` arguments on modal/static/transient
entry points are legacy z-axis spin values in revolutions/s and
revolutions/s^2. The Timoshenko element path converts them internally to rad/s
and rad/s^2 before adding them to the element rigid-body angular state. Use
`rbData` when a full angular velocity or angular acceleration vector is needed.

## Loads

Structural loads may be nodal, distributed per unit length, follower loads, or
generalized reduced-order loads. Validation tests should name the load category
because the expected reaction sign and magnitude differ.

External point loads passed as `Fexternal` and `Fdof` use global DOF numbering.
Concentrated nodal terms parsed by [`applyConcentratedTerms`](@ref) are stored
on the model and included during element assembly. Distributed element loads are
mapped through element transforms and shape functions.

## Gravity And Spin Terms

`FEAModel` controls gravity, spin-up, nonlinear, and aeroelastic options.
Spin-softening, centrifugal, and Coriolis effects depend on the selected analysis
path. A structural-only test should state whether gravity and rotor speed are
enabled.

The current public steady `staticAnalysis` scalar `Omega` and `OmegaDot` inputs
are z-axis spin terms supplied in revolutions/s and revolutions/s^2. The
Timoshenko element path converts them internally to rad/s and rad/s^2 before
adding them to the element rigid-body angular state. Global x-axis spin for HAWT
rotor-only motion still needs an explicit public frame/state owner before it can
replace the scalar z-axis path.

`el.rotationalEffects[i] != 1` disables scalar and vector spin terms for element
`i`. Use this for non-rotating tower or support elements in a mixed mesh.

## Solver Flow

Most workflows follow the same assembly path:

1. [`FEAModel`](@ref) constructs joint transforms, reduced DOF lists, boundary
   condition maps, nonlinear parameters, and default nodal terms.
2. [`initialElementCalculations`](@ref) performs the first Timoshenko element
   pass and caches element storage for reuse.
3. `TimoshenkoMatrixWrap!` loops over elements, builds element inputs, calls the
   Timoshenko element, and assembles global matrices/vectors.
4. Joint constraints are applied with `jointTransform' * K * jointTransform` and
   `jointTransform' * F` when constrained joints are present.
5. Boundary conditions are applied by row/column modification for static and
   transient solves, or by row/column removal for modal solves.
6. The selected solver path computes displacements, eigenvalues, or reduced
   modal coordinates.
7. Strains and reaction forces are recovered after the solve.

## Solver Paths

| Path | Entry point | Main behavior |
| --- | --- | --- |
| Modal | [`modal`](@ref) | Assemble `K`, `M`, `C`, reduce for joints/BCs, build state-space matrix, solve eigenvalues, normalize mode shapes. |
| Static | [`staticAnalysis`](@ref) | Apply load stepping, solve linear or nonlinear displacement updates, recover strain and reactions. |
| Transient | [`structuralDynamicsTransient`](@ref) | Build an effective Newmark-beta or Dean system, iterate nonlinear displacements when enabled, update displacement/velocity/acceleration. |
| ROM | [`reducedOrderModel`](@ref), [`structuralDynamicsTransientROM`](@ref) | Project full matrices into modal coordinates and time-integrate reduced states. |

For nonlinear direct iteration (`"DI"`), the solver replaces the displacement
state and optionally relaxes the update. For Newton-Raphson (`"NR"`), the solver
uses displacement increments where that tangent path is implemented. Some
low-level nonlinear selective-stiffness paths intentionally reject unsupported
Newton-Raphson calls; those guardrails are pinned in `test/HelperFunctions.jl`.
