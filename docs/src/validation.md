# Validation and Testing

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
eccentric beams under scalar z-axis spin acceleration. The eccentric case also
documents the current root-moment gap: the x-force from the offset is present,
while the reported spin-axis torque omits the corresponding `y^2` contribution.

## Acceptance Rules

- Pin exact matrix entries for utility functions.
- Use physically named displacement, strain, reaction, and frequency channels in
  solver tests.
- State whether a test is a regression, a unit/API contract, or a validation
  against another solver.
- Keep generated modal matrices and output files ignored or written to a
  temporary path.

## Current Known Gaps

- `steady.jl`, `intermediate.jl`, and parts of `timoshenko.jl` still have lower
  coverage than the utility layer.
- The transient GXBeam comparison reports large RMS differences for some
  channels; those are regression diagnostics, not proof of close physical
  agreement.
- Several public-looking functions are unexported and need clearer docstrings
  before the API should be considered stable.
