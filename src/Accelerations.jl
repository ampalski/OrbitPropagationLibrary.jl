# Will need to pass options, if not the whole input, through p
function force_model(x, p, t)
    accel = SA[0.0, 0, 0]
    jd0 = p[1]
    opts = p[2]
    jd = JDate(SA[jd0.epoch[1], jd0.epoch[2] + t / 86400.0], jd0.system)
    r_sun = SA[0.0, 0, 0]
    r_sc_sun = SA[0.0, 0, 0]

    if opts.third_body_sun || opts.solar_radiation_pressure
        r_sun = _sun_pos(jd)
        r_sc_sun = r_sun - @view x[pos_inds]
    end

    if opts.third_body_moon
        r_moon = _moon_pos(jd)
        r_sc_moon = r_moon - @view x[pos_inds]
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
        @views a_srp = _srpaccel(x[pos_inds], r_sun, opts)
        accel += a_srp
    end
    if opts.non_spherical
        P = SMatrix{3, 3}(mod2j200076_matrix(jd))
        N = SMatrix{3, 3}(tod2mod76_matrix(jd))
        R = SMatrix{3, 3}(pef2tod76_matrix(jd))
        W = SMatrix{3, 3}(itrf2pef76_matrix(jd))
        pos_ecef = W' * R' * N' * P' * @view x[1:3]
        a_nonsph = _nonsph_accel(pos_ecef, opts.degree, opts.order)
        accel += P * N * R * W * a_nonsph
    end

    r = norm(@view x[pos_inds])
    a_twobody = -μ / r^3 * @view x[pos_inds]
    accel += a_twobody

    return [x[vel_inds];accel]
end

#from Montenbruck & Gill's "Satellite Orbits"
function _sun_pos(JD)
    JDTDB = convert_jd(JD, :TDB)
    AS2RAD = 2.0 * pi / 360 / 3600
    epsilon = 23.43929111 * pi / 180.0     # Obliquity of J2000 ecliptic
    T = ((JDTDB.epoch[1] - 2451545.0) + JDTDB.epoch[2]) / 36525

    temp = 282.94 * pi / 180
    M = 357.5256 + 35999.049 * T
    M *= pi / 180
    λs = (6892 * sin(M) + 72 * sin(2 * M)) * AS2RAD
    λs += temp + M
    rs = 1.0e6 * (149.619 - 2.499 * cos(M) - 0.021 * cos(2 * M))

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
function _moon_pos(JD)
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


function _srpaccel(r, r_sun, opts)
    # Check for umbra & penumbra conditions
    r_sc_sun = r_sun - r
    shadow_val = _shadowfraction(r, r_sun)
    # If umbra, return zeros
    if shadow_val == 0.0
        return zeros(3)
    end
    p_srp = 1367.0 / 3.0e8 # W*s/m3
    a_srp = -p_srp * opts.coefficient_of_radiation *
        opts.area / opts.mass * AU^2 / norm(r_sc_sun)^3 * r_sc_sun
    a_srp /= 1000 # convert to km
    return shadow_val * a_srp
    # If penumbra, return partial acceleration value
end

function _factorial_term(l, m)
    δk = m == 0 ? 1 : 2
    temp = δk * (2 * l + 1)
    for i in (l - m + 1):(l + m)
        temp /= i
    end
    return sqrt(temp)
    #return factorial(l - m) * δk * (2 * l + 1) / factorial((l + m)))
end

function _nonsph_accel(pos::AbstractVector, degree::Int, order::Int)
    lmax = degree
    mmax = order
    r = norm(pos)
    V = spzeros(lmax + 3, mmax + 3)
    W = spzeros(lmax + 3, mmax + 3)
    V[1, 1] = REarth / r
    W[1, 1] = 0.0

    # Zonals
    Rr2 = V[1, 1] / r
    V[2, 1] = Rr2 * pos[3] * REarth / r
    for l in 2:(lmax + 1)
        V[l + 1, 1] = (2 * l - 1) / (l) * pos[3] * Rr2 * V[l, 1]
        V[l + 1, 1] -= (l - 1) / l * Rr2 * REarth * V[l - 1, 1]
    end

    # Tesserals and Sectorials
    for m in 1:(mmax + 1)
        V[m + 1, m + 1] = (2 * m - 1) * (pos[1] * Rr2 * V[m, m] - pos[2] * Rr2 * W[m, m])
        W[m + 1, m + 1] = (2 * m - 1) * (pos[1] * Rr2 * W[m, m] + pos[2] * Rr2 * V[m, m])
        for l in (m + 1):(lmax + 1)
            V[l + 1, m + 1] = (2 * l - 1) / (l - m) * pos[3] * Rr2 * V[l, m + 1]
            V[l + 1, m + 1] -= (l + m - 1) / (l - m) * Rr2 * REarth * V[l - 1, m + 1]

            W[l + 1, m + 1] = (2 * l - 1) / (l - m) * pos[3] * Rr2 * W[l, m + 1]
            W[l + 1, m + 1] -= (l + m - 1) / (l - m) * Rr2 * REarth * W[l - 1, m + 1]
        end
    end

    #Calculate accelerations
    ax = 0.0
    ay = 0.0
    az = 0.0

    for m in mmax:-1:0 # backwards to try to fix precision stuff
        for l in lmax:-1:m
            # Grab and un-normalize the coefficients
            norm = _factorial_term(l, m)
            C = NormGravityModel_C[l + 1, m + 1] * norm
            S = NormGravityModel_S[l + 1, m + 1] * norm

            if m == 0
                ax -= C * V[l + 2, 2]
                ay -= C * W[l + 2, 2]
                az -= (l + 1) * C * V[l + 2, 1]
            else
                term = (l - m + 1) * (l - m + 2)
                ax += 0.5 * (
                    (-C * V[l + 2, m + 2] - S * W[l + 2, m + 2]) +
                        term * (C * V[l + 2, m] + S * W[l + 2, m])
                )
                ay += 0.5 * (
                    (-C * W[l + 2, m + 2] + S * V[l + 2, m + 2]) +
                        term * (-C * W[l + 2, m] + S * V[l + 2, m])
                )
                az += (l - m + 1) * (-C * V[l + 2, m + 1] - S * W[l + 2, m + 1])
            end
        end
    end
    return μ / REarth^2 * SA[ax, ay, az]
end
