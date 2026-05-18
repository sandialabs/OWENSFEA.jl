# Model Assembly

This page documents the contracts that are easiest to break during refactors.

## Mesh And Elements

`Mesh` stores node coordinates, element connectivity, element types, mesh
segments, and optional structural blade/tower indexing. `El` stores element-level
properties and distributed loads. Keep node ordering and element orientation
stable because those choices define local beam frames and reaction signs.

## Section Properties

`SectionPropsArray` carries spanwise structural properties:

- `EA`, `GJ`, `EIyy`, `EIzz`, and optional coupling terms;
- mass per unit length and inertias;
- center-of-mass and aerodynamic-center offsets;
- twist and aerodynamic/structural reference offsets.

When section data comes from `OWENSPreComp`, document whether offsets are
relative to the leading edge, reference axis, shear center, tension center, or
center of mass before mapping into FEA inputs.

## Boundary Conditions

Boundary conditions are expressed as node number, local DOF number, and
prescribed value. Utility tests pin the reduced DOF map and static boundary
condition application. A constrained DOF should either be removed from the solve
or assigned the prescribed value consistently in the full displacement vector.

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

## Concentrated Terms

`applyConcentratedTerms` parses concentrated mass, stiffness, damping, and load
terms into global arrays. Four-column rows apply diagonal terms, while
five-column rows specify a row/column DOF pair. Loads are one-dimensional and
ignore the second DOF in five-column input.

Use temporary files for parser tests and prefer direct `data = ...` input for
unit tests that only need to pin matrix assembly.
