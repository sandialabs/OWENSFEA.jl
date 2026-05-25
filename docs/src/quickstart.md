# Quickstart

OWENSFEA workflows are easiest to debug when mesh, sections, boundary
conditions, and loads are assembled explicitly. The cantilever tests provide the
smallest maintained examples:

- `test/CantileverBeamModal.jl` for modal analysis;
- `test/CantileverBeamDisplacement.jl` for nonlinear/static displacement;
- `test/UnsteadyBeam.jl` for transient dynamics and GXBeam comparison;
- `test/HelperFunctions.jl` for joint, boundary-condition, and nodal-term
  utility contracts.

## Install

```julia
using Pkg
Pkg.add(PackageSpec(url = "https://github.com/sandialabs/OWENSFEA.jl.git"))
```

For local toolkit development:

```julia
using Pkg
Pkg.develop(path = ".")
Pkg.instantiate()
```

The docs project uses local package sources and should be built with Julia 1.11
or newer.

## Minimal Structural Workflow

1. Build or read a `Mesh`.
2. Build an `El` object containing element-level material and section data.
3. Create a `FEAModel` with analysis type, boundary conditions, joint data, and
   solver options.
4. Run the desired analysis path: modal, steady, unsteady, or ROM.
5. Check displacement, strain, reaction, and modal outputs in the expected frame.

The package exports selected utilities, but the high-level model types are often
used through qualified names:

```julia
import OWENSFEA

fea = OWENSFEA.FEAModel(; analysisType = "M", numNodes = 2)
```

## First Checks

Before trusting a result:

- confirm the mesh has six structural degrees of freedom per node;
- verify constrained DOF maps and joint transforms;
- verify section properties are in SI units;
- check whether rotations are radians internally or degrees in a file;
- compare reactions and displacement signs against a hand calculation for a
  single-load case.
