```@meta
CurrentModule = OWENSFEA
```

# API Map

This page maps the main user-facing API. Complete docstrings are grouped by
source file on [Autodocs by Source](@ref).

## Model Data

| Name | Purpose |
| --- | --- |
| [`FEAModel`](@ref) | Solver options, boundary conditions, joint transforms, nonlinear parameters, and nodal terms. |
| `Mesh` | Node coordinates, connectivity, element type metadata, and structural indexing. |
| [`Ort`](@ref) | Element orientation helper data. |
| [`SectionPropsArray`](@ref) | Spanwise section stiffness, mass, offsets, aerodynamic helper fields, and shear stiffness inputs. |
| [`El`](@ref) | Per-element section properties, lengths, angles, and rotational-effects flags. |
| [`DispData`](@ref) | Displacement, velocity, acceleration, and previous-step state input. |
| [`DispOut`](@ref) | Displacement, velocity, acceleration, and optional modal state output. |
| [`ElStrain`](@ref) | Element strain and curvature channels. |
| [`ROM`](@ref) | Reduced-order matrices and coefficient storage. |

## Analysis Entry Points

| Name | Purpose |
| --- | --- |
| [`modal`](@ref) | Linearized modal/state-space analysis. |
| [`autoCampbellDiagram`](@ref) | RPM sweep for OWENSFEA or GXBeam modal frequencies. |
| [`staticAnalysis`](@ref) | Static/steady displacement, strain, and reaction solve. |
| [`initialElementCalculations`](@ref) | Cache initial element data before repeated analyses. |
| [`structuralDynamicsTransient`](@ref) | Full-order transient structural dynamics. |
| [`reducedOrderModel`](@ref) | Build reduced-order model matrices. |
| [`structuralDynamicsTransientROM`](@ref) | Transient structural dynamics in reduced coordinates. |

## Assembly And Mapping

| Name | Purpose |
| --- | --- |
| [`createJointTransform`](@ref) | Build reduced/full DOF transform from joint rows. |
| [`calculateReducedDOFVector`](@ref) | Build a reduced DOF list from constraints. |
| [`constructReducedDispVectorMap`](@ref) | Build reduced displacement maps with BCs. |
| [`calculateBCMap`](@ref) | Map prescribed BCs into equation numbering. |
| [`applyBC`](@ref) | Apply static/transient prescribed boundary conditions. |
| [`assembly!`](@ref) | Assemble element matrix/vector into global arrays. |
| [`assemblyMatrixOnly`](@ref) | Assemble an element matrix without a load vector. |
| [`applyConcentratedTerms`](@ref) | Parse concentrated mass, stiffness, damping, and load terms. |

## Element And Postprocessing Helpers

| Name | Purpose |
| --- | --- |
| [`defaultTransverseShearStiffness`](@ref) | Compute default Timoshenko transverse shear stiffness. |
| [`transverseShearStiffness`](@ref) | Interpolate explicit or default shear stiffness at a quadrature point. |
| [`calculateTimoshenkoElementInitialRun`](@ref) | Initial Timoshenko element calculation. |
| [`calculateTimoshenkoElementNL`](@ref) | Nonlinear/current-state Timoshenko element calculation. |
| [`calculateTimoshenkoElementStrain`](@ref) | Recover element strain channels. |
| [`calculateTimoshenkoElementNLSS`](@ref) | Nonlinear selective-stiffness helper pinned by utility tests. |
| [`calculateStructureMassProps`](@ref) | Integrate stored element mass properties. |
| [`calculateReactionForceAtNode`](@ref) | Recover reaction load at a node. |
| [`calculateStrainForElements`](@ref) | Recover strains for all elements. |
| [`findElementsAssociatedWithNodeNumber`](@ref) | Locate element associations for a node. |
| [`getGP`](@ref) | Gauss quadrature points and weights. |
| [`calculateShapeFunctions`](@ref) | Linear/quadratic Lagrange shape functions and Jacobian. |

## Load Stepping And Modal Helpers

| Name | Purpose |
| --- | --- |
| [`updateLoadStep`](@ref) | Update static load-step state after an iteration. |
| [`adaptiveLoadStepping`](@ref) | Adaptive load-step branch logic. |
| [`applyBCModal`](@ref) | Remove constrained rows/columns for modal matrices. |
| [`extractFreqDamp`](@ref) | Convert eigenpairs to frequency, damping, and normalized mode-shape data. |
| [`constructReducedDispVecFromEigVec`](@ref) | Reconstruct reduced displacement vector from an eigenvector and BC map. |

## Index

```@index
```
