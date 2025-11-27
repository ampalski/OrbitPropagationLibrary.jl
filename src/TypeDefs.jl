const validFrames = [:ITRF, :PEF, :TOD, :TEME, :MOD, :J2000]

export OplIn, InitialState, OplOut
abstract type OplIn end
abstract type InitialState end
abstract type OplOut end

export Cart_InitialState
struct Cart_InitialState <: InitialState
    state::SVector{6,Float64}
    epoch::JulianDate
    frame::Symbol
end

export build_cartesian_state
"""
    state0 = build_cartesian_state(stateVec, epoch, frame)

Builds an InitialState object for the `stateVec` specified at `epoch` in `frame`

Allowable frames are: `:ITRF`, `:PEF`, `:TOD`, `:TEME`, `:MOD`, `:J2000`

The input value for `epoch` is the Julian Date returned in two pieces, in the 
usual SOFA manner, which is designed to preserve time resolution. The full 
Julian Date is available as a single number by adding the two components of the
vector.

# Inputs
* `state::AbstractVector`: 6-element vector cartesian state vector (position and velocity)
* `epoch::JulianDate`: The Julian Date for the given state
* `frame::Symbol`: The frame the state is specified in

"""
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
    Δt::Float64
    outputs::Vector{Symbol}
    frame::Symbol
    finalTime::JulianDate
end

export build_base_output
"""
    output = build_base_output(Δt, outputs, frame, finalTime)

Builds an OPLOut object for the given options.

The `outputs` argument is a vector of the specific data points desired at each
output timestep, specified as Symbols. Allowable options are currently: 
* Cartesian positions: :x, :y, :z, :pos
* Cartesian velocities: :vx, :vy, :vz, :vel
* Cartesian state: :state
* Classical orbital elements: :coe, :sma, :ecc, :inc, :Ω, :ω, :M, :ν (nu)

Output timesteps are given every `Δt` from the initial state until `finalTime`.
A value of `Δt=0` will use the default propagator timesteps.
A value of `Δt=-1` will only output data at `finalTime`

Allowable frames are: `:ITRF`, `:PEF`, `:TOD`, `:TEME`, `:MOD`, `:J2000`

The input value for `epoch` is the Julian Date returned in two pieces, in the 
usual SOFA manner, which is designed to preserve time resolution. The full 
Julian Date is available as a single number by adding the two components of the
vector.

# Inputs
* `Δt::Real`: Specify the timestep between output values
* `outputs::Vector{Symbol}`: Specify the output data points.
* `frame::Symbol`: Specify the frame the `outputs` are in.
* `finalTime::JulianDate`: Specify the final time to propagate to.

"""
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
    return BaseOplOut(Float64(Δt), outputs, frame, finalTime)
end

export TwoBody_OplIn
struct TwoBody_OplIn <: OplIn
    state0::InitialState
    output::OplOut
end

export build_twobody_input
"""
    input = build_twobody(state0, output)

Builds an OplIn object for the given settings.

The three inputs are best created using the helper functions:
* `state0` using `build_cartesian_state`
* `output` using `build_base_output`

# Inputs
* `state0::InitialState`: Object describing time and state information at the beginning of propagation
* `output::OplOut`: Object describing the desired output, including final propagation time and what data to return.

"""
function build_twobody_input(state0::InitialState, output::OplOut)
    return TwoBody_OplIn(state0, output)
end

export NumericalOptions
struct NumericalOptions
    non_spherical::Bool
    third_body_moon::Bool
    third_body_sun::Bool
    thrust::Bool
    solar_radiation_pressure::Bool
    drag::Bool
    area::Float64
    mass::Float64
    coefficient_of_radiation::Float64
    degree::Int
    order::Int
end

# in the documentation for this function, need to include units
export build_numerical_options
"""
    options = build_numerical_options(; kwargs...)

Builds a NumericalOptions object for the given key word options.

# Inputs
* `use_non_spherical::Bool=false`: If true, utilize terms up to `degree` and `order`
* `use_third_body_sun::Bool=false`: Third body effects from the sun.
* `use_third_body_moon::Bool=false`: Third body effects from the moon.
* `use_solar_radiation_pressure=false`: If true, utilize `area`, `mass`, and `coefficient_of_radiation`
* `use_drag::Bool=false`: Not currently implemented.
* `use_thrust::Bool=false`: Not currently implemented.
* `area::Float64=1.0`: Average area presented to solar flux, in m²
* `mass:Float64=1000.0`: Mass of the spacecraft, in kg
* `coefficient_of_radiation::Float64=1.0`: C_R term
* `degree::Int=20`: Maximum degree of the geopotential model, up to 200
* `order::Int=20`: Maximum order of the geopotential model, up to 200

"""
function build_numerical_options(;
    use_non_spherical::Bool=false,
    use_third_body_sun::Bool=false,
    use_third_body_moon::Bool=false,
    use_solar_radiation_pressure::Bool=false,
    use_drag::Bool=false,
    use_thrust::Bool=false,
    area::Float64=1.0,
    mass::Float64=1000.0,
    coefficient_of_radiation::Float64=1.0,
    degree::Int=20,
    order::Int=20,
)
    # For any options turned on, check that the corresponding settings are
    # in place

    if area < 0
        error("Negative spacecraft areas are not possible")
    end
    if mass < 0
        error("Negative spacecraft masses are not possible")
    end
    if degree < 0 || order < 0
        error("Geopotential model degree and order must be non-negative")
    end
    if degree > 200 || order > 200
        error("Geopotential model degree and order must be ≤ 200")
    end

    # Build the struct
    return NumericalOptions(
        use_non_spherical, use_third_body_moon,
        use_third_body_sun, use_thrust,
        use_solar_radiation_pressure, use_drag,
        area, mass, coefficient_of_radiation, degree, order)
end

export Numerical_OplIn
struct Numerical_OplIn <: OplIn
    state0::InitialState
    options::NumericalOptions
    output::OplOut
end

export build_numerical_input
"""
    input = build_numerical_input(state0, options, output)

Builds an OplIn object for the given settings.

The three inputs are best created using the helper functions:
* `state0` using `build_cartesian_state`
* `options` using `build_numerical_options`
* `output` using `build_base_output`

# Inputs
* `state0::InitialState`: Object describing time and state information at the beginning of propagation
* `options::NumericalOptions`: Object describing the force model and spacecraft settings to use for propagation
* `output::OplOut`: Object describing the desired output, including final propagation time and what data to return.

"""
function build_numerical_input(
    state0::InitialState,
    options::NumericalOptions,
    output::OplOut,
)
    return Numerical_OplIn(state0, options, output)
end

