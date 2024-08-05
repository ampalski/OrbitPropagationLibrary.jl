# TODO: accept a state vector, it's frame, and the output request and frame
# Covert to the desired elements, stash into a DataFrame

function constructoutput(
    rf::AbstractVector,
    vf::AbstractVector,
    output::OplOut,
)
    # Initialize data frame
    df = DataFrame()
    for col in output.outputs
        df[!, col] = _getoutputtype(col)
    end

    # Build each data type
    #TODO: loop over time first, add row by row
    row = [_convertoutput(rf, vf, col) for col in output.outputs]
    push!(df, row)
    return df
end

function _getoutputtype(type::Symbol)
    if type == :x || type == :y || type == :z
        return Float64[]
    elseif type == :vx || type == :vy || type == :vz
        return Float64[]
    elseif type == :state || type == :coe
        return Vector{SVector{6,Float64}}()
    elseif type == :pos || type == :vel
        return Vector{SVector{3,Float64}}()
    elseif type == :sma || type == :ecc || type == :inc || type == :Ω
        return Float64[]
    elseif type == :ω || type == :M || type == :ν
        return Float64[]
    end
end

function _convertoutput(rf, vf, type)
    if type == :pos
        return SA[rf...]
    elseif type == :vel
        return SA[vf...]
    elseif type == :state
        return SA[vcat(rf, vf)...]
    end
    #TODO: Add the rest
end

# const validOutputs = [:x, :y, :z, :vx, :vy, :vz, :state, :pos, :vel, :coe,
#     :sma, :ecc, :inc, :Ω, :ω, :M, :ν]
