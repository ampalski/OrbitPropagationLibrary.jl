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
    elseif type == :x
        return rf[1]
    elseif type == :y
        return rf[2]
    elseif type == :z
        return rf[3]
    elseif type == :vx
        return vf[1]
    elseif type == :vy
        return vf[2]
    elseif type == :vz
        return vf[3]
    elseif type == :coe || type == :sma || type == :ecc || type == :inc ||
           type == :Ω || type == :ω || type == :M || type == :ν
        coes = state_to_classical(rf, vf)
        if type == :coe
            return coes
        elseif type == :sma
            return coes[1]
        elseif type == :ecc
            return coes[2]
        elseif type == :inc
            return coes[3]
        elseif type == :Ω
            return coes[4]
        elseif type == :ω
            return coes[5]
        elseif type == :M
            return coes[6]
        elseif type == :ν
            return coes[7]
        end
    end

end

# const validOutputs = [:x, :y, :z, :vx, :vy, :vz, :state, :pos, :vel, :coe,
#     :sma, :ecc, :inc, :Ω, :ω, :M, :ν]
