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

    dt = 86400.0 * ((tf.epoch[1] - t0.epoch[1]) + (tf.epoch[2] - t0.epoch[2]))
    if dt < 0
        error("Output time is before input time")
    end

    # Convert input to whatever Universal expects
    # TODO: convert frames as needed
    # TODO: use the dt in the output, and loop over to get intermediate points
    rf, vf = universalkepler(input.state0.state[1:3],
        input.state0.state[4:6], dt, μ)

    # Construct output
    return constructoutput(rf, vf, input.output)

    # return [rf; vf]
end
