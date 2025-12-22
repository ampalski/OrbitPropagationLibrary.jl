# include specific stuff here
include("UnivKepler.jl")
include("Accelerations.jl")
include("Outputs.jl")

export propagate
"""
    output_df = propagate(input)

Propagate the satellite state specified in `input` and provide an output 
dataframe.

The `input` struct is best created through the helper functions 
`build_twobody_input` or `build_numerical_input`. Currently only Two-Body and
Numerical are supported.

# Inputs
* `input::OplIn`: Input struct specifying input, options, and desired output.

# Examples
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
"""
function propagate(input::OplIn)
    return if typeof(input) != TwoBody_OplIn || typeof(input) != Numerical_OplIn
        error("Unsupported input type")
    end
end

function propagate(input::TwoBody_OplIn)
    # make sure output time is past input time
    tf = input.output.finalTime
    if tf.system != input.state0.epoch.system
        t0 = convert_jd(input.state0.epoch, tf.system)
    else
        t0 = input.state0.epoch
    end

    totalTime = 86400.0 * ((tf.epoch[1] - t0.epoch[1]) + (tf.epoch[2] - t0.epoch[2]))
    if totalTime < 0
        error("Output time is before input time")
    end

    # Convert input to whatever Universal expects
    if input.output.Δt > 0
        dt = input.output.Δt
        T = collect(dt:dt:totalTime)
        if T[end] != totalTime
            push!(T, totalTime)
        end
    else
        T = [totalTime]
    end
    rf = Vector{SVector{3, Float64}}()
    vf = Vector{SVector{3, Float64}}()
    r = SA[input.state0.state[1:3]...]
    v = SA[input.state0.state[4:6]...]
    if input.state0.frame != :J2000
        converted_state = convert_state([r; v], input.state0.frame, :J2000, t0)
        r = SA[converted_state[1:3]...]
        v = SA[converted_state[4:6]...]
    end
    for i in eachindex(T)
        dt = i == 1 ? T[1] : T[i] - T[i - 1]

        r, v = universalkepler(r, v, dt, μ)
        push!(rf, r)
        push!(vf, v)
    end

    # Construct output
    return _constructoutput(rf, vf, T, input)

    # return [rf; vf]
end

function propagate(input::Numerical_OplIn)
    # make sure output time is past input time
    tf = input.output.finalTime
    if tf.system != input.state0.epoch.system
        t0 = convert_jd(input.state0.epoch, tf.system)
    else
        t0 = input.state0.epoch
    end

    totalTime = 86400.0 * ((tf.epoch[1] - t0.epoch[1]) + (tf.epoch[2] - t0.epoch[2]))
    if totalTime < 0
        error("Output time is before input time")
    end

    if input.output.Δt > 0
        dt = input.output.Δt
        T = collect(dt:dt:totalTime)
        if T[end] != totalTime
            push!(T, totalTime)
        end
    else
        T = [totalTime]
    end
    rf = Vector{SVector{3, Float64}}()
    vf = Vector{SVector{3, Float64}}()
    r = SA[input.state0.state[1:3]...]
    v = SA[input.state0.state[4:6]...]
    if input.state0.frame != :J2000
        converted_state = convert_state([r; v], input.state0.frame, :J2000, t0)
        r = SA[converted_state[1:3]...]
        v = SA[converted_state[4:6]...]
    end
    # for i in eachindex(T)
    #     dt = i == 1 ? T[1] : T[i] - T[i-1]
    #
    #     r, v = universalkepler(r, v, dt, μ)
    #     push!(rf, r)
    #     push!(vf, v)
    # end
    # should be able to use ODEProblem to directly get the times I want
    ode_t = (0, totalTime)
    opts = (input.state0.epoch, input.options)
    prob = ODEProblem(force_model, [r; v], ode_t, opts)
    # sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12, saveat=dt)
    sol = solve(prob, Tsit5(), reltol = 1.0e-12, abstol = 1.0e-12).(T)
    #need to test the saveat version against using the solution as a function
    for x in sol
        push!(rf, x[1:3])
        push!(vf, x[4:6])
    end

    # Construct output
    return _constructoutput(rf, vf, T, input)

    # return [rf; vf]
end
