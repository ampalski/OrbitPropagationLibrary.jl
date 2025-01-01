# include specific stuff here
include("UnivKepler.jl")
include("Accelerations.jl")
include("Outputs.jl")

export propagate
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
    rf = Vector{SVector{3,Float64}}()
    vf = Vector{SVector{3,Float64}}()
    r = SA[input.state0.state[1:3]...]
    v = SA[input.state0.state[4:6]...]
    if input.state0.frame != :J2000
        converted_state = convert_state([r; v], input.state0.frame, :J2000, t0)
        r = SA[converted_state[1:3]...]
        v = SA[converted_state[4:6]...]
    end
    for i in eachindex(T)
        dt = i == 1 ? T[1] : T[i] - T[i-1]

        r, v = universalkepler(r, v, dt, μ)
        push!(rf, r)
        push!(vf, v)
    end

    # Construct output
    return constructoutput(rf, vf, T, input)

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
    rf = Vector{SVector{3,Float64}}()
    vf = Vector{SVector{3,Float64}}()
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
    opts = [input.state0.epoch, input.options]
    prob = ODEProblem(force_model, [r; v], ode_t, opts)
    # sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12, saveat=dt)
    sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12).(T)
    #need to test the saveat version against using the solution as a function
    for x in sol
        push!(rf, x[1:3])
        push!(vf, x[4:6])
    end

    # Construct output
    return constructoutput(rf, vf, T, input)

    # return [rf; vf]
end
