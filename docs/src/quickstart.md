# Quick Start

```@meta
CurrentModule = OWENSFEA
```

OWENSFEA workflows are easiest to debug when mesh, sections, boundary
conditions, and loads are assembled explicitly. The cantilever tests provide the
smallest maintained examples:

- `test/CantileverBeamModal.jl` for modal analysis;
- `test/CantileverBeamDisplacement.jl` for nonlinear/static displacement;
- `test/UnsteadyBeam.jl` for transient dynamics and GXBeam comparison;
- `test/HelperFunctions.jl` for joint, boundary-condition, and nodal-term
  utility contracts.

## Install Or Develop

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

## Minimal Cantilever Model

The following is the smallest direct-constructor pattern used by the tests. It
builds a straight cantilever with six degrees of freedom per node, fixed root
boundary conditions, and constant Timoshenko section properties.

```julia
import OWENSFEA
using LinearAlgebra

L = 0.5
nelem = 8
nnodes = nelem + 1
num_dof_per_node = 6

x = collect(range(0.0, L; length = nnodes))
y = zeros(nnodes)
z = zeros(nnodes)
conn = hcat(collect(1:nelem), collect(2:nnodes))

mesh = OWENSFEA.Mesh(
    collect(1:nnodes),
    nelem,
    nnodes,
    x,
    y,
    z,
    collect(1:nelem),
    conn,
    zeros(Int, nelem),
    [nelem],
    zeros(1, 1),
    zeros(Int, 1, 1),
    zeros(Int, 1, 1),
)

width = 0.05
height = 0.02
area = width * height
E = 2.1e11
nu = 0.28
G = E / (2 * (1 + nu))
density = 7800.0
Iyy = width * height^3 / 12
Izz = width^3 * height / 12
J = Iyy + Izz

function constant_section()
    zero2 = [0.0, 0.0]
    GA = fill(5 / 6 * G * area, 2)
    return OWENSFEA.SectionPropsArray(
        zero2,
        zero2,
        fill(density * area, 2),
        fill(E * Iyy, 2),
        fill(E * Izz, 2),
        fill(G * J, 2),
        fill(E * area, 2),
        fill(density * Iyy, 2),
        fill(density * Izz, 2),
        fill(density * J, 2),
        zero2,
        zero2,
        zero2,
        zero2,
        zero2,
        zero2,
        zero2,
        zero2,
        zero2,
        zero2,
        zero2,
        zero2,
        zero2,
        zero2,
        nothing,
        nothing,
        zero2,
        zero2,
        GA,
        GA,
    )
end

section_props = [constant_section() for _ in 1:nelem]
el = OWENSFEA.El(
    section_props,
    fill(L / nelem, nelem),
    zeros(nelem),
    zeros(nelem),
    zeros(nelem),
    ones(nelem),
)

pBC = [
    1 1 0
    1 2 0
    1 3 0
    1 4 0
    1 5 0
    1 6 0
]
joint = zeros(0, 8)
```

The high-level model types are normally used through qualified names:

```julia
feamodel = OWENSFEA.FEAModel(;
    analysisType = "M",
    dataOutputFilename = "none",
    joint,
    pBC,
    numNodes = mesh.numNodes,
    gravityOn = false,
    numModes = 8,
)
```

## Modal Analysis

```julia
modal_result = OWENSFEA.modal(feamodel, mesh, el)
freq = modal_result[1]
damp = modal_result[2]
first_mode_hz = freq[1]
```

`test/CantileverBeamModal.jl` compares this style of model against GXBeam and
Euler-Bernoulli bending frequencies. The test extracts every other OWENSFEA
state-space mode when comparing to GXBeam because the eigenvalues appear in
conjugate pairs.

## Static Tip Load

```julia
feamodel = OWENSFEA.FEAModel(;
    analysisType = "TNB",
    dataOutputFilename = "none",
    joint,
    pBC,
    numNodes = mesh.numNodes,
    nlOn = false,
    gravityOn = false,
    iterationType = "LINEAR",
)

el_storage = OWENSFEA.initialElementCalculations(feamodel, el, mesh)
displ0 = zeros(mesh.numNodes * num_dof_per_node)
tip_z_dof = (mesh.numNodes - 1) * num_dof_per_node + 3

displ, strain, success, reaction = OWENSFEA.staticAnalysis(
    feamodel,
    mesh,
    el,
    displ0,
    0.0,
    0.0,
    el_storage;
    Fdof = [tip_z_dof],
    Fexternal = [100.0],
)
```

`staticAnalysis` returns the full displacement vector, element strain objects, a
success flag, and reaction loads for each node/element slot used by the current
implementation. For a cantilever sanity check, the linear tip deflection should
scale with `P * L^3 / (3 * E * I)`.

## Transient Step

```julia
disp_data = OWENSFEA.DispData(displ0, zero(displ0), zero(displ0), zero(displ0))
dt = 0.002
time = dt

strain, disp_out, reaction = OWENSFEA.structuralDynamicsTransient(
    feamodel,
    mesh,
    el,
    disp_data,
    0.0,
    0.0,
    time,
    dt,
    el_storage,
    [100.0],
    [tip_z_dof],
    I(3),
    zeros(9),
)
```

For a time march, update [`DispData`](@ref) from `disp_out.displ_sp1`,
`disp_out.displdot_sp1`, and `disp_out.displddot_sp1` before the next step.

## First Checks

Before trusting a result:

- confirm the mesh has six structural degrees of freedom per node;
- verify constrained DOF maps and joint transforms;
- verify section properties are in SI units;
- check whether rotations are radians internally or degrees at an IO boundary;
- compare reactions and displacement signs against a hand calculation for a
  single-load case.
