```@meta
CurrentModule = OWENSFEA
```

# Autodocs by Source

This page lists package docstrings. Prefer the [API Map](@ref) for a
user-facing overview.

## Source Groups

- `structs.jl`: model data containers and constructors.
- `modal.jl`: modal analysis, Campbell diagrams, and modal helpers.
- `steady.jl`: static analysis and load stepping.
- `unsteady.jl`: transient full-order dynamics and matrix mapping helpers.
- `rom.jl`: reduced-order model generation and integration.
- `timoshenko.jl`: element-level Timoshenko calculations.
- `utilities.jl`: joints, BCs, assembly, concentrated terms, reactions, strains,
  and shape functions.
- `intermediate.jl`: analysis-type dispatch into element assembly.

## Complete Autodocs

```@autodocs
Modules = [OWENSFEA]
Order = [:type, :function]
```
