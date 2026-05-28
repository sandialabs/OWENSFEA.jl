# Validation and Testing

```@meta
CurrentModule = OWENSFEA
```

OWENSFEA has three useful test layers:

| Layer | Test files | What is pinned |
| --- | --- | --- |
| Helper/API | `test/HelperFunctions.jl` | quadrature, shape functions, reduced DOF maps, joint transforms, concentrated terms, and boundary-condition application. |
| Static/modal beams | `test/CantileverBeam*.jl` | beam displacement, modal frequencies, rotations, and nonlinear/static behavior. |
| Transient comparison | `test/UnsteadyBeam.jl` | time response compared with GXBeam for the maintained cantilever case. |

The helper layer also pins an off-axis nonlinear selective-stiffness matrix from
`calculateTimoshenkoElementNLSS`. That case is a regression guard for future
strain-stiffening and nonlinear ROM work: it fixes the current direct-iteration
matrix entries and requires unsupported Newton-Raphson calls to fail fast with a
clear error instead of reaching undefined tangent-matrix state.

`test/SpinningBeamDiagnostics.jl` pins structural-only spin reactions before
gyric or spin-softening equation changes. The cases include straight and
eccentric beams under scalar z-axis spin acceleration, sectional-CG offsets
under steady scalar z-axis spin, plus a steady rigid-body spin-vector path for
non-z spin axes. The sectional-CG case verifies that rotating-frame inertial
loads are evaluated at the center of mass so the offset force and offset moment
cancel about the spin axis for steady spin.

## Maintained Physical Checks

| Check | Source | Expected behavior |
| --- | --- | --- |
| Straight cantilever modes | `test/CantileverBeamModal.jl` | OWENSFEA matches GXBeam modes within 1 percent and Euler-Bernoulli bending modes within 5 percent for the documented modes. |
| Static tip deflection | `test/CantileverBeamDisplacement.jl` | Linear response follows `P * L^3 / (3 * E * Iyy)` and top-fiber bending strain follows the analytical beam relation. |
| Swept rotating beam | `test/CantileverBeamRotating*.jl` | 45 degree swept-beam response and rotating modal trends stay close to GXBeam for the tested RPM range. |
| Transient tip force | `test/UnsteadyBeam.jl` | Newmark-beta transient root force/moment channels stay inside current broad regression tolerances against GXBeam. |
| Nonlinear ROM | `test/NonlinearROM.jl` | ROM tip history remains within 5 percent of the full nonlinear transient result for the maintained cantilever. |
| Spin reactions | `test/SpinningBeamDiagnostics.jl` | Steady spin, spin acceleration, and center-of-mass offset reactions match hand-derived force and moment balances for pinned cases. |

## Running Tests

From the package root:

```julia
using Pkg
Pkg.test()
```

Or from a shell:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Some tests compare against GXBeam and can be slower than utility tests. When
debugging a small utility change, run the specific included test file first,
then run the full suite before merging.

## Acceptance Rules

- Pin exact matrix entries for utility functions.
- Use physically named displacement, strain, reaction, and frequency channels in
  solver tests.
- State whether a test is a regression, a unit/API contract, or a validation
  against another solver.
- Keep generated modal matrices and output files ignored or written to a
  temporary path.
- Include units, frame, load direction, and enabled solver options in any new
  physics test.

## Current Known Gaps

- `steady.jl`, `intermediate.jl`, and parts of `timoshenko.jl` still have lower
  coverage than the utility layer.
- The transient GXBeam comparison reports large RMS differences for some
  channels; those are regression diagnostics, not proof of close physical
  agreement.
- Several public-looking functions are unexported and need clearer docstrings
  before the API should be considered stable.
- There is no separate `examples/` directory in this package at present; test
  files are the maintained examples.
