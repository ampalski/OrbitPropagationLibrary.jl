function wrapto2pi(input::Float64)
    return rem2pi(input, RoundDown)
end

function wraptopi(input::Float64)
    return rem2pi(input, RoundNearest)
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

function unit(a::AbstractVector)
    return a ./ norm(a)
end

function anglevec(r1::AbstractVector, r2::AbstractVector)
    return acos(dot(r1, r2) / norm(r1) / norm(r2))
end

#Requires radians
function Rx(angle::Real)
    c = cos(angle)
    s = sin(angle)

    return SMatrix{3, 3, Float64}(
        [
            +1.0 0.0 0.0;
            0.0 +c +s;
            0.0 -s +c
        ]
    )
end

function _shadowfraction(r, r_sun)
    α_umb = 0.004609793064071904 #rad
    α_pen = 0.004695061438837353

    # Same side of earth as sun, no shadow
    if r' * r_sun > 0
        return 1.0
    end
    # simplified x-y coordinates
    ξ = anglevec(r, -r_sun)
    rn = norm(r)
    horiz = rn * cos(ξ)
    vert = rn * sin(ξ)
    # outer edge of the penumbra region
    x = REarth / sin(α_pen)
    pen_vert = tan(α_pen) * (x + horiz)
    # check if outside penumbra region
    if vert > pen_vert
        return 1.0
    end
    # outer edge of the umbra cone
    y = REarth / sin(α_umb)
    umb_vert = tan(α_pen) * (y - horiz)
    # check if inside umbra cone
    if vert < umb_vert
        return 0.0
    end
    # fractional shadow
    # return (vert - umb_vert) / (pen_vert - umb_vert)
    # Below version calculates occulting discs, from Montenbruck
    a = asin(RSun / norm(r_sun - r))
    b = asin(REarth / rn)
    c = acos((-r' * r_sun - r) / (rn * norm(r_sun - r)))
    x = (c^2 + a^2 - b^2) / (2 * c)
    y = sqrt(a^2 - x^2)
    A = a^2 * acos(x / a) + b^2 * acos((c - x) / b) - c * y
    return 1 - A / (π * a^2)
    #TODO: Worth checking the simple version against the occulting discs,
    # and verify that the conditions on Montenbruck pg 83 hold true.
end

function kronecker(a::Real, b::Real)
    if a == b
        return 1.0
    else
        return 0.0
    end
end
