@testset "TwoBody" begin
    jd0, _ = datevec2jdate([2020.0, 1, 10, 19, 0, 0])
    jdf, _ = datevec2jdate([2020.0, 1, 11, 19, 0, 0])

    state0 = build_cartesian_state(
        [33333.0, 22222, -44444, 0.40731, 2.191364, 1.401164],
        jd0,
        :J2000,
    )
    output = build_base_output(-1.0, [:state], :J2000, jdf)
    input1 = build_twobody_input(state0, output)
    out1 = propagate(input1).state[1]

    options = build_numerical_options()
    input2 = build_numerical_input(state0, options, output)
    out2 = propagate(input2).state[1]

    # Make sure analytic and numerical match
    @test isapprox(out1[1], out2[1], atol=0.0001)
    @test isapprox(out1[2], out2[2], atol=0.0001)
    @test isapprox(out1[3], out2[3], atol=0.0001)
    @test isapprox(out1[4], out2[4], atol=0.0001)
    @test isapprox(out1[5], out2[5], atol=0.0001)
    @test isapprox(out1[6], out2[6], atol=0.0001)

    # Make sure they match the STK output
    @test isapprox(-37349.726871, out2[1], atol=0.1)
    @test isapprox(-40380.983615, out2[2], atol=0.1)
    @test isapprox(34121.543544, out2[3], atol=0.1)
    @test isapprox(0.108347, out2[4], atol=0.001)
    @test isapprox(-1.596218, out2[5], atol=0.001)
    @test isapprox(-1.834136, out2[6], atol=0.001)
end

@testset "SRP" begin
    jd0, _ = datevec2jdate([2020.0, 1, 10, 19, 0, 0])
    jdf, _ = datevec2jdate([2020.0, 1, 11, 19, 0, 0])

    state0 = build_cartesian_state(
        [33333.0, 22222, -44444, 0.40731, 2.191364, 1.401164],
        jd0,
        :J2000,
    )
    output = build_base_output(-1.0, [:state], :J2000, jdf)
    options = build_numerical_options(
        use_solar_radiation_pressure=true,
        area=1.0,
        mass=50.0,
    )
    input = build_numerical_input(state0, options, output)
    out = propagate(input).state[1]

    # Make sure the results match the STK output
    @test isapprox(-37349.956682, out[1], atol=0.1)
    @test isapprox(-40380.688951, out[2], atol=0.1)
    @test isapprox(34122.160974, out[3], atol=0.1)
    @test isapprox(0.108333, out[4], atol=0.001)
    @test isapprox(-1.596226, out[5], atol=0.001)
    @test isapprox(-1.834115, out[6], atol=0.001)
end

@testset "sun" begin
    jd0, _ = datevec2jdate([2020.0, 1, 10, 19, 0, 0])
    jdf, _ = datevec2jdate([2020.0, 1, 11, 19, 0, 0])

    state0 = build_cartesian_state(
        [33333.0, 22222, -44444, 0.40731, 2.191364, 1.401164],
        jd0,
        :J2000,
    )
    output = build_base_output(-1.0, [:state], :J2000, jdf)
    options = build_numerical_options(
        use_third_body_sun=true,
    )
    input = build_numerical_input(state0, options, output)
    out = propagate(input).state[1]

    # Make sure the results match the STK output
    @test isapprox(-37356.460243, out[1], atol=0.1)
    @test isapprox(-40376.046875, out[2], atol=0.1)
    @test isapprox(34131.799469, out[3], atol=0.2) #?
    @test isapprox(0.108135, out[4], atol=0.001)
    @test isapprox(-1.596294, out[5], atol=0.001)
    @test isapprox(-1.833754, out[6], atol=0.001)
end

@testset "moon" begin
    jd0, _ = datevec2jdate([2020.0, 1, 10, 19, 0, 0])
    jdf, _ = datevec2jdate([2020.0, 1, 11, 19, 0, 0])

    state0 = build_cartesian_state(
        [33333.0, 22222, -44444, 0.40731, 2.191364, 1.401164],
        jd0,
        :J2000,
    )
    output = build_base_output(-1.0, [:state], :J2000, jdf)
    options = build_numerical_options(
        use_third_body_moon=true,
    )
    input = build_numerical_input(state0, options, output)
    out = propagate(input).state[1]

    # Make sure the results match the STK output
    @test isapprox(-37366.834068, out[1], atol=0.1)
    @test isapprox(-40373.381366, out[2], atol=0.1)
    @test isapprox(34137.534728, out[3], atol=0.1)
    @test isapprox(0.107948, out[4], atol=0.001)
    @test isapprox(-1.596264, out[5], atol=0.001)
    @test isapprox(-1.833407, out[6], atol=0.001)
end

@testset "NonSpherical" begin
    jd0, _ = datevec2jdate([2020.0, 1, 10, 19, 0, 0])
    jdf, _ = datevec2jdate([2020.0, 1, 11, 19, 0, 0])

    state0 = build_cartesian_state(
        [33333.0, 22222, -44444, 0.40731, 2.191364, 1.401164],
        jd0,
        :J2000,
    )
    output = build_base_output(-1.0, [:state], :J2000, jdf)
    options = build_numerical_options(
        use_non_spherical=true,
        degree=20,
        order=20,
    )
    input = build_numerical_input(state0, options, output)
    out = propagate(input).state[1]

    # Make sure the results match the STK output
    @test isapprox(-37351.471751, out[1], atol=0.1)
    @test isapprox(-40379.503769, out[2], atol=0.1)
    @test isapprox(34123.924202, out[3], atol=0.1)
    @test isapprox(0.108247, out[4], atol=0.001)
    @test isapprox(-1.596257, out[5], atol=0.001)
    @test isapprox(-1.834056, out[6], atol=0.001)
end

@testset "FullUp" begin
    jd0, _ = datevec2jdate([2020.0, 1, 10, 19, 0, 0])
    jdf, _ = datevec2jdate([2020.0, 1, 11, 19, 0, 0])

    state0 = build_cartesian_state(
        [33333.0, 22222, -44444, 0.40731, 2.191364, 1.401164],
        jd0,
        :J2000,
    )
    output = build_base_output(-1.0, [:state], :J2000, jdf)
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
    input = build_numerical_input(state0, options, output)
    out = propagate(input).state[1]

    # Make sure the results match the STK output
    @test isapprox(-37375.542895, out[1], atol=0.1)
    @test isapprox(-40366.665574, out[2], atol=0.1)
    @test isapprox(34150.786305, out[3], atol=0.2)
    @test isapprox(0.107621, out[4], atol=0.001)
    @test isapprox(-1.596387, out[5], atol=0.001)
    @test isapprox(-1.832925, out[6], atol=0.001)
end
