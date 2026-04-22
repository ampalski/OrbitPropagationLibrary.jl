# OrbitPropagationLibrary

[![Build Status](https://github.com/ampalski/OrbitPropagationLibrary.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ampalski/OrbitPropagationLibrary.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

An Earth-centered propagation library built for realistic accelerations, helpful interfaces, and flexibility for mission analysis and mission planning.

This package is at minimum viable product, allowing for two-body propagation as well as a useful set of additional accelerations, to include non-spherical Earth effects, third-body effects from the sun and moon, and solar radiation pressure. Additional accelerations (drag, tides, etc.), maneuver handling, and more will be built in to future releases.

In order to maintain compatibility with the USSF standards, initial implementation focuses on IAU-76 and other pre-2000 models, with the exception of the EGM-2006 model.

## Usage Example

```
using OrbitPropagationLibrarySOFA
using OrbitPropagationLibrary
jd0, _ = datevec2jdate([2020.0, 1, 10, 19, 0, 0])
jdf, _ = datevec2jdate([2020.0, 1, 11, 19, 0, 0])

state0 = build_cartesian_state(
    [33333.0, 22222, -44444, 0.40731, 2.191364, 1.401164],
    jd0,
    :J2000,
)
output_config = build_base_output(-1.0, [:state], :J2000, jdf)
input1 = build_twobody_input(state0, output_config)
out1 = propagate(input1)

options = build_numerical_options(
    use_non_spherical=true,
    use_third_body_sun=true,
    use_third_body_moon=true,
    use_solar_radiation_pressure=true,
    area=1.0,
    mass=50.0,
    degree=20,
    order=20,
)
input2 = build_numerical_input(state0, options, output_config)
out2 = propagate(input2)
```

## Future Plans
- Additional accelerations, to include relativity and tides.
- Maneuver handling.
- Interfacing for mission optimization.

