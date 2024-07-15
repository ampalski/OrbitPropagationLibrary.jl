const validFrames = [:ITRF, :PEF, :TOD, :TEME, :MOD, :J2000]

export OplIn, InitialState, OplOut
abstract type OplIn end
abstract type InitialState end
abstract type OplOut end

export Cart_InitialState
struct Cart_InitialState <: InitialState
    cartState::SVector{6,Float64}
    epoch::JulianDate
    frame::Symbol
end

export build_cartesian_state
function build_cartesian_state(
    state::AbstractVector,
    epoch::JulianDate,
    frame::Symbol,
)
    if !(frame in validFrames)
        error("Invalid coordinate frame provided.")
    end
    if state isa SVector{6,Float64}
        return Cart_InitialState(state, epoch, frame)
    end

    if length(state) != 6
        error("State vectors must be of length 6, with stacked position and velocity vectors")
    end

    return Cart_InitialState(SA[state...], epoch, frame)
end

const validOutputs = [:x, :y, :z, :vx, :vy, :vz, :state, :pos, :vel, :coe,
    :sma, :ecc, :inc, :Ω, :ω, :M, :ν]

export BaseOplOut
# Δt governs how often output is produced. 0 uses any default values from the
# propagator (if any), -1 produces only the final values, other positive values
# give the timestep in seconds
struct BaseOplOut <: OplOut
    Δt::Real
    outputs::Vector{Symbol}
    frame::Symbol
    finalTime::JulianDate
end

export build_base_output
function build_base_output(
    Δt::Real,
    outputs::Vector{Symbol},
    frame::Symbol,
    finalTime::JulianDate,
)
    if !(frame in validFrames)
        error("Invalid coordinate frame provided.")
    end
    for output in outputs
        if !(output in validOutputs)
            error("$(output) is not a valid output selection.")
        end
    end
    return BaseOplOut(Δt, outputs, frame, finalTime)
end

export TwoBody_OplIn
struct TwoBody_OplIn <: OplIn
    state0::InitialState
    output::OplOut
end

export build_twobody_input
function build_twobody_input(state0::InitialState, output::OplOut)
    return TwoBody_OplIn(state0, output)
end


