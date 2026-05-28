# Model Assembly

```@meta
CurrentModule = OWENSFEA
```

This page documents the contracts that are easiest to break during refactors.

## Required Objects

Most workflows use the same four objects:

| Object | Role |
| --- | --- |
| `Mesh` | Node coordinates, element connectivity, element types, segment metadata, and optional rotating/non-rotating metadata. |
| [`SectionPropsArray`](@ref) | Two-point spanwise section properties for each element. |
| [`El`](@ref) | Per-element section property, length, orientation, twist/roll, and rotational-effects flag arrays. |
| [`FEAModel`](@ref) | Analysis type, boundary conditions, joints, nonlinear parameters, damping, gravity/spin flags, and nodal terms. |

Create `Mesh`, `SectionPropsArray`, and `El` first, then construct
[`FEAModel`](@ref) with `numNodes = mesh.numNodes`. The model constructor builds
joint transforms, reduced-DOF maps, boundary-condition maps, and empty nodal
term containers when those are not supplied.

## Mesh And Elements

`Mesh` stores node coordinates, element connectivity, element types, mesh
segments, and optional structural blade/tower indexing. `El` stores element-level
properties and distributed loads. Keep node ordering and element orientation
stable because those choices define local beam frames and reaction signs.

The maintained tests build straight and swept cantilever meshes directly in
`test/testdeps.jl`. That helper computes element `Psi_d` and `Theta_d` from the
node-to-node vector, stores lengths in `Ort.Length`, and uses joint rows to weld
the two beam segments together.

Element DOFs are node-major in public vectors:

```text
node i: ux, uy, uz, theta_x, theta_y, theta_z
```

For a node `n` and local DOF `d`, the global DOF is `(n - 1) * 6 + d`.

## Section Properties

`SectionPropsArray` carries spanwise structural properties:

- `EA`, `GJ`, `EIyy`, `EIzz`, optional transverse shear stiffnesses
  `GAy`/`GAz`, and optional coupling terms;
- mass per unit length and inertias;
- center-of-mass and aerodynamic-center offsets;
- twist and aerodynamic/structural reference offsets.

When section data comes from `OWENSPreComp`, document whether offsets are
relative to the leading edge, reference axis, shear center, tension center, or
center of mass before mapping into FEA inputs.

When `GAy` and `GAz` are omitted, the Timoshenko element computes both
transverse shear stiffnesses from `EA`, Poisson's ratio 0.3, and a 5/6 shear
correction unless `poisson_ratio` or `shear_correction` arrays are supplied on
the section properties. Composite or externally generated section-property
workflows should pass explicit `GAy` and `GAz` values when available; otherwise
they should pass the material Poisson ratio and shear-correction factor used to
derive the transverse shear stiffnesses.

The tests use two common section-property patterns:

- rectangular-beam properties with `rhoA`, `EA`, `EIyy`, `EIzz`, `GJ`,
  `rhoIyy`, `rhoIzz`, and `rhoJ` set from beam dimensions;
- explicit `GAy` and `GAz` for rotating modal comparisons, so OWENSFEA and
  GXBeam use the same transverse shear stiffness.

Most properties are passed as two values per element endpoint and interpolated
at quadrature points. Use SI units unless the caller's file parser clearly
documents a conversion.

## Boundary Conditions

Boundary conditions are expressed as node number, local DOF number, and
prescribed value. Utility tests pin the reduced DOF map and static boundary
condition application. A constrained DOF should either be removed from the solve
or assigned the prescribed value consistently in the full displacement vector.

Root-fixed cantilever tests use:

```julia
pBC = [
    1 1 0
    1 2 0
    1 3 0
    1 4 0
    1 5 0
    1 6 0
]
```

`makeBCdata` builds [`BC_struct`](@ref) maps from this matrix, the reduced DOF
list, and any joint transform. Static solves apply prescribed values directly
with [`applyBC`](@ref); modal solves remove constrained rows/columns.

## Joints

Joint rows use this contract:

| Column | Meaning |
| --- | --- |
| 1 | joint number |
| 2 | master node |
| 3 | slave node |
| 4 | joint type |
| 5 | joint mass or rigid-bar `lx`, depending on context |
| 6 | rigid-bar `ly` |
| 7 | rigid-bar `lz` or `psi`, depending on context |
| 8 | `theta` |

Joint types currently covered by unit tests:

- `0`: fixed/welded;
- `1`: pinned translational constraint;
- `2`: hinge along local 2 axis;
- `3`: hinge along local 1 axis;
- `4`: hinge along local 3 axis;
- `5`: rigid-bar constraint.

The rigid-bar path couples slave translations to master rotations through the
bar offset. Tests pin this transform because sign mistakes directly affect
reaction forces.

If a model has no joints, pass `zeros(0, 8)` or another zero-row matrix rather
than a placeholder welded joint. The constructor accepts legacy placeholders,
but explicit empty joint data is less ambiguous.

## Concentrated Terms

`applyConcentratedTerms` parses concentrated mass, stiffness, damping, and load
terms into global arrays. Four-column rows apply diagonal terms, while
five-column rows specify a row/column DOF pair. Loads are one-dimensional and
ignore the second DOF in five-column input.

Use temporary files for parser tests and prefer direct `data = ...` input for
unit tests that only need to pin matrix assembly.

Supported data tags are:

| Tag | Meaning |
| --- | --- |
| `"M"` | diagonal concentrated mass |
| `"K"` | diagonal concentrated stiffness |
| `"C"` | diagonal concentrated damping |
| `"F"` | concentrated nodal load |
| `"M6"` | mass matrix entry with row/column DOFs |
| `"K6"` | stiffness matrix entry with row/column DOFs |
| `"C6"` | damping matrix entry with row/column DOFs |

Example direct input from the helper-function tests:

```julia
nodalinputdata = Any[
    1 "M" 2 1.5
    1 "K" 3 2.5
    1 "C" 4 3.5
    1 "F" 5 4.5
]

nodal_terms = OWENSFEA.applyConcentratedTerms(1, 6; data = nodalinputdata)
```

## Analysis Types

`FEAModel.analysisType` selects solver-side behavior:

| Value | Use |
| --- | --- |
| `"M"` | modal analysis and linearized matrices |
| `"S"` | internal static/steady path set by [`staticAnalysis`](@ref) |
| `"TNB"` | Newmark-beta transient dynamics |
| `"TD"` | Dean transient dynamics |
| `"ROM"` | reduced-order transient dynamics |
| `"stiff"` | pass-through structural path used by coupled workflows that do not solve FEA |
| `"GX"` | GXBeam branch in [`autoCampbellDiagram`](@ref), not an OWENSFEA element solve |

Use the entry-point function as the authority when possible. For example,
[`staticAnalysis`](@ref) sets the model analysis type to `"S"` internally, and
[`structuralDynamicsTransientROM`](@ref) sets it to `"ROM"`.
