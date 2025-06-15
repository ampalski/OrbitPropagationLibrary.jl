function wrapto2pi(input::Float64)
    p2 = 2 * pi
    while input < 0
        input += p2
    end

    return input % p2
end

function wraptopi(input::Float64)
    p2 = 2 * pi
    w = input % (p2)
    if abs(w) > 0.5
        w -= sign(input) * p2
    end
    return w
end

function cross(a::AbstractVector, b::AbstractVector)
    if !(length(a) == length(b) == 3)
        throw(DimensionMismatch("cross product is only defined for vectors of length 3"))
    end
    a1, a2, a3 = a
    b1, b2, b3 = b

    return [a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1]
end

function dot(a::AbstractVector, b::AbstractVector)
    if length(a) != length(b)
        throw(DimensionMismatch("dot product requires input vectors of the same length"))
    end

    return a' * b
end

function norm(a::AbstractVector)
    return sqrt(dot(a, a))
end

function anglevec(r1::AbstractVector, r2::AbstractVector)
    return acos(dot(r1, r2) / norm(r1) / norm(r2))
end

#Requires radians
function Rx(angle::Real)


    c = cos(angle)
    s = sin(angle)

    return [+1.0 0.0 0.0;
        0.0 +c +s;
        0.0 -s +c]
end
