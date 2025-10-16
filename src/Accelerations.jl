# Will need to pass options, if not the whole input, through p
function force_model(x, p, t)
    accel = zeros(3)
    jd0 = p[1]
    opts = p[2]
    jd = JDate(SA[jd0.epoch[1], jd0.epoch[2]+t/86400.0], jd0.system)

    if opts.third_body_sun || opts.solar_radiation_pressure
        r_sun = sun_pos(jd)
        r_sc_sun = r_sun - x[1:3]
    end

    if opts.third_body_moon
        r_moon = moon_pos(jd)
        r_sc_moon = r_moon - x[1:3]
        temp1 = norm(r_sc_moon)^3
        temp2 = norm(r_moon)^3
        a_moon = μMOON * (r_sc_moon ./ temp1 - r_moon ./ temp2)
        accel += a_moon
    end
    if opts.third_body_sun
        temp1 = norm(r_sc_sun)^3
        temp2 = norm(r_sun)^3
        a_sun = μSUN * (r_sc_sun ./ temp1 - r_sun ./ temp2)
        accel += a_sun
    end
    if opts.solar_radiation_pressure
        a_srp = srpaccel(x[1:3], r_sun, opts)
        accel += a_srp
    end

    r = norm(x[1:3])
    a_twobody = -μ / r^3 * x[1:3]
    accel += a_twobody

    return SA[x[4:6]...; accel...]
end


#from Montenbruck & Gill's "Satellite Orbits"
export sun_pos
function sun_pos(JD)
    JDTDB = convert_jd(JD, :TDB)
    AS2RAD = 2.0 * pi / 360 / 3600
    epsilon = 23.43929111 * pi / 180.0     # Obliquity of J2000 ecliptic
    T = ((JDTDB.epoch[1] - 2451545.0) + JDTDB.epoch[2]) / 36525

    temp = 282.94 * pi / 180
    M = 357.5256 + 35999.049 * T
    M *= pi / 180
    λs = (6892 * sin(M) + 72 * sin(2 * M)) * AS2RAD
    λs += temp + M
    rs = 1e6 * (149.619 - 2.499 * cos(M) - 0.021 * cos(2 * M))

    pos = rs * Rx(-epsilon) * [cos(λs), sin(λs), 0.0]

    return pos
end

# From vallado, deprecated in favor of the above
# function sun_pos(JD)
#     JDUT1 = convert_jd(JD, :UT1)
#     JDTDB = convert_jd(JD, :TDB)
#     T = ((JDUT1.epoch[1] - 2451545.0) + JDUT1.epoch[2]) / 36525
#     TTDB = ((JDTDB.epoch[1] - 2451545.0) + JDTDB.epoch[2]) / 36525
#
#     ecliptic_longitude = 280.46 + 36000.771 * T
#     MSUN = 357.5277233 + 35999.05034 * TTDB
#     ecliptic_longitude += (1.914666471 * sind(MSUN) +
#                            0.019994643 * sind(2 * MSUN))
#
#     rsun = 1.000140612 - 0.016708617 * cosd(MSUN) -
#            0.000139589 * cosd(2 * MSUN)
#     rsun *= 149597870.7
#     ecliptic_obl = 23.439291 - 0.01461 * TTDB
#
#     r1 = rsun * cosd(ecliptic_longitude)
#     r2 = rsun * cosd(ecliptic_obl) * sind(ecliptic_longitude)
#     r3 = rsun * sind(ecliptic_obl) * sind(ecliptic_longitude)
#
#     return [r1, r2, r3]
# end

# From Vallado, not as accurate to JPL Ephems as the below version
# function moon_pos(JD)
#     JDTDB = convert_jd(JD, :TDB)
#     T = ((JDTDB.epoch[1] - 2451545.0) + JDTDB.epoch[2]) / 36525
#     Rearth = 6378.1363
#     ecliptic_longitude = 218.32 + 481267.8813 * T +
#                          6.29 * sind(134.9 + 477198.85 * T) -
#                          1.27 * sind(259.2 - 413335.38 * T) +
#                          0.66 * sind(235.7 + 890534.23 * T) +
#                          0.21 * sind(269.9 + 954397.7 * T) -
#                          0.19 * sind(357.5 + 35999.05 * T) -
#                          0.11 * sind(186.6 + 966404.05 * T)
#     ecliptic_latitude = 5.13 * sind(93.3 + 483202.03 * T) +
#                         0.28 * sind(228.2 + 960400.87 * T) -
#                         0.28 * sind(318.3 + 6003.18 * T) -
#                         0.17 * sind(217.6 - 407332.2 * T)
#     horizontal_parallax = 0.9508 + 0.0518 * cosd(134.9 + 477198.85 * T) +
#                           0.0095 * cosd(259.2 - 413335.38 * T) +
#                           0.0078 * cosd(235.7 + 890534.23 * T) +
#                           0.0028 * cosd(269.9 + 954397.7 * T)
#     ecliptic_obl = 23.439291 - 0.0130042 * T - 1.64e-7 * T^2 + 5.04e-7 * T^3
#     rmoon = Rearth / sind(horizontal_parallax)
#     r1 = rmoon * cosd(ecliptic_latitude) * cosd(ecliptic_longitude)
#     r2 = rmoon * (cosd(ecliptic_obl) * cosd(ecliptic_latitude) *
#                   sind(ecliptic_longitude) - sind(ecliptic_obl) *
#                                              sind(ecliptic_latitude))
#     r3 = rmoon * (sind(ecliptic_obl) * cosd(ecliptic_latitude) *
#                   sind(ecliptic_longitude) + cosd(ecliptic_obl) *
#                                              sind(ecliptic_latitude))
#     return [r1, r2, r3]
# end

# From Montenbruck & Gill's "Satellite Orbits"
function moon_pos(JD)
    JDTDB = convert_jd(JD, :TDB)
    T = ((JDTDB.epoch[1] - 2451545.0) + JDTDB.epoch[2]) / 36525
    epsilon = 23.43929111 * pi / 180.0     # Obliquity of J2000 ecliptic
    AS2RAD = 2.0 * pi / 360 / 3600

    # Mean elements of lunar orbit, units in radians
    L0 = pi / 180 * (218.31617 + 481267.88088 * T − 1.3972 * T)
    l = pi / 180 * (134.96292 + 477198.86753 * T)
    lp = pi / 180 * (357.52543 + 35999.04944 * T)
    F = pi / 180 * (93.27283 + 483202.01873 * T)
    D = pi / 180 * (297.85027 + 445267.11135 * T)

    # Moon longitude
    λm = 22640.0 * sin(l) + 769 * sin(2 * l) -
         4586 * sin(l - 2 * D) + 2370 * sin(2 * D) -
         668 * sin(lp) - 412 * sin(2 * F) -
         212 * sin(2 * l - 2 * D) - 206 * sin(l + lp - 2 * D) +
         192 * sin(l + 2 * D) - 165 * sin(lp - 2D) +
         148 * sin(l - lp) - 125 * sin(D) -
         110 * sin(l + lp) - 55 * sin(2 * F - 2 * D)
    λm *= AS2RAD
    λm += L0

    # Moon latitude
    temp = 412 * sin(2 * F) + 541 * sin(lp)
    temp *= AS2RAD
    βm = 18520 * sin(F + λm - L0 + temp) -
         526 * sin(F - 2 * D) + 44 * sin(l + F - 2 * D) -
         31 * sin(F - l - 2D) - 25 * sin(F - 2 * l) -
         23 * sin(lp + F - 2 * D) + 21 * sin(F - l) +
         11 * sin(F - lp - 2 * D)
    βm *= AS2RAD

    # Moon distance (in km)
    rm = 385000.0 - 20905 * cos(l) - 3699 * cos(2 * D - l) -
         2956 * cos(2 * D) - 570 * cos(2 * l) + 246 * cos(2 * l - 2 * D) -
         205 * cos(lp - 2 * D) - 171 * cos(l + 2 * D) -
         152 * cos(l + lp - 2 * D)

    # Position
    pos = rm * Rx(-epsilon) * [cos(λm) * cos(βm), sin(λm) * cos(βm), sin(βm)]

    return pos
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

export srpaccel
function srpaccel(r, r_sun, opts)
    # Check for umbra & penumbra conditions
    r_sc_sun = r_sun - r
    shadow_val = _shadowfraction(r, r_sun)
    # If umbra, return zeros
    if shadow_val == 0.0
        return zeros(3)
    end
    p_srp = 1367.0 / 3e8 # W*s/m3 
    a_srp = -p_srp * opts.coefficient_of_radiation *
            opts.area / opts.mass * AU^2 / norm(r_sc_sun)^3 * r_sc_sun
    a_srp /= 1000 # convert to km
    # TODO: doesn't match results on pg 606
    return shadow_val * a_srp
    # If penumbra, return partial acceleration value
end
