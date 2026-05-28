# Developer Guide

```@meta
CurrentModule = OWENSFEA
```

This guide is for changes to OWENSFEA internals, tests, and documentation.

## Work From The Tests

Read the relevant test before changing solver logic:

| Area | Start with |
| --- | --- |
| shape functions, assembly, BCs, joints, nodal terms | `test/HelperFunctions.jl` |
| modal frequencies and matrix export | `test/CantileverBeamModal.jl` |
| static/nonlinear displacement and strain | `test/CantileverBeamDisplacement.jl` |
| spin and rotating-frame reactions | `test/SpinningBeamDiagnostics.jl` |
| swept rotating beams and Campbell trends | `test/CantileverBeamRotating*.jl` |
| transient dynamics | `test/UnsteadyBeam.jl` |
| reduced-order dynamics | `test/NonlinearROM.jl` |

The tests are also the maintained example source. When adding docs, prefer a
small excerpt or distilled pattern from a test over an unverified new example.

## Numerical Contracts

Treat these as first-class implementation constraints:

- public displacement and force vectors use node-major six-DOF ordering;
- `pBC` values are prescribed displacements/rotations, not penalty terms;
- joint constraints are enforced through transform matrices;
- section properties use SI units and two endpoint values per element;
- mesh orientation fields are in degrees, while internal angular states are in
  radians;
- scalar `Omega`/`OmegaDot` are legacy z-axis rev/s inputs, while `rbData`
  angular entries are rad/s vectors;
- `rotationalEffects` controls whether an element receives spin terms;
- generated matrices and output files should go to temporary paths in tests.

Prefer explicit equations, transforms, and constraints over tuning penalties.
If a change alters units, signs, nondimensionalization, load scaling, or state
updates, add a small test that exposes the convention.

## Solver Changes

Before editing `steady.jl`, `unsteady.jl`, `intermediate.jl`, `rom.jl`, or
`timoshenko.jl`, identify which path owns the behavior:

| Path | Core files |
| --- | --- |
| element matrices and strains | `src/timoshenko.jl` |
| assembly dispatch by analysis type | `src/intermediate.jl` |
| static load stepping | `src/steady.jl` |
| transient Newmark-beta/Dean integration | `src/unsteady.jl` |
| modal state-space eigensolve | `src/modal.jl` |
| reduced-order projection/integration | `src/rom.jl` |
| maps, BCs, joints, concentrated terms, reactions | `src/utilities.jl` |

Keep low-level utility tests exact and high-level physics tests tolerant enough
to allow harmless numerical reordering. For modal changes, compare mode pairing
and ordering carefully because state-space eigenvalues are paired.

## Adding Public API

Many useful functions are currently unexported. Before exporting another name:

- add or update its docstring;
- add it to the reference map if it is user-facing;
- add a focused test for argument shape, units, and frame expectations;
- avoid mutating [`FEAModel`](@ref) fields as a hidden side effect unless the
  existing entry point already does so and the behavior is documented.

If a function is only meant for package internals, keep it unexported and document
it through the grouped autodocs rather than promoting it in examples.

## Documentation Workflow

Build docs from the package root:

```bash
julia --project=docs docs/make.jl
```

The docs project points `OWENSFEA` at the local package via `docs/Project.toml`,
so the build should not require network access when dependencies are already
instantiated.

When adding examples:

- use plain `julia` code blocks for longer examples that should not run during
  docs build;
- use `@example` only for short snippets that are cheap and deterministic;
- include the source test file near every test-backed example;
- avoid documenting a result unless it is actually asserted in tests or produced
  by a checked script.

## Review Checklist

Before considering a solver or docs change complete:

- run the narrow test or docs build that exercises the change;
- inspect numerical outputs for units and sign consistency;
- update captions/text when a result changes;
- keep generated artifacts out of version control unless they are intentionally
  tracked;
- check `git status` and avoid reverting unrelated edits from other work.
