# include specific stuff here
include("UnivKepler.jl")
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
    # TODO: convert frames as needed
    if input.output.Δt > 0
        dt = input.output.Δt
        T = collect(dt:dt:totalTime)
        if T[end] != totalTime
            push!(T, totalTime)
        end
    else
        T = [totalTime]
    end
    rf = Vector{Vector{Float64}}()
    vf = Vector{Vector{Float64}}()
    r = input.state0.state[1:3]
    v = input.state0.state[4:6]
    for i in eachindex(T)
        dt = i == 1 ? T[1] : T[i] - T[i-1]

        r, v = universalkepler(r, v, dt, μ)
        push!(rf, r)
        push!(vf, v)
    end

    # Construct output
    return constructoutput(rf, vf, T, input.output)

    # return [rf; vf]
end
