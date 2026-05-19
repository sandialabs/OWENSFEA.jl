# Theory, Frames, and Units

## Beam Coordinates

OWENSFEA uses six degrees of freedom per node: three translations and three
rotations. Element-level stiffness, mass, and load terms are assembled in local
beam coordinates and mapped into the global system through element orientation
data.

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
| bending stiffness | N m^2 |
| torsional stiffness | N m^2 |
| axial stiffness | N |

## Rotations

Internal rotations and angular velocities are in radians. Some mesh inputs,
orientation helpers, and plots use degrees for readability; conversion should
occur at the IO or plotting boundary.

## Loads

Structural loads may be nodal, distributed per unit length, follower loads, or
generalized reduced-order loads. Validation tests should name the load category
because the expected reaction sign and magnitude differ.

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
