# Will need to pass options, if not the whole input, through p
function force_model(x, p, t)
    accel = zeros(3)
    jd0 = p[1]
    opts = p[2]
    jd = JDate(SA[jd0.epoch[1], jd0.epoch[2]+t/86400.0], jd0.system)

    if opts.third_body_moon
        r_moon = moon_pos(jd)
        r_sc_moon = r_moon - x[1:3]
        temp1 = norm(r_sc_moon)^3
        temp2 = norm(r_moon)^3
        a_moon = μMOON * (r_sc_moon ./ temp1 - r_moon ./ temp2)
        accel += a_moon
    end
    if opts.third_body_sun
        r_sun = sun_pos(jd)
        r_sc_sun = r_sun - x[1:3]
        temp1 = norm(r_sc_sun)^3
        temp2 = norm(r_sun)^3
        a_sun = μSUN * (r_sc_sun ./ temp1 - r_sun ./ temp2)
        accel += a_sun
    end

    r = norm(x[1:3])
    a_twobody = -μ / r^3 * x[1:3]
    accel += a_twobody

    return SA[x[4:6]...; accel...]
end

function moon_pos(JD)
    JDTDB = convert_jd(JD, :TDB)
    T = ((JDTDB.epoch[1] - 2451545.0) + JDTDB.epoch[2]) / 36525
    Rearth = 6378.1363
    ecliptic_longitude = 218.32 + 481267.8813 * T +
                         6.29 * sind(134.9 + 477198.85 * T) -
                         1.27 * sind(259.2 - 413335.38 * T) +
                         0.66 * sind(235.7 + 890534.23 * T) +
                         0.21 * sind(269.9 + 954397.7 * T) -
                         0.19 * sind(357.5 + 35999.05 * T) -
                         0.11 * sind(186.6 + 966404.05 * T)
    ecliptic_latitude = 5.13 * sind(93.3 + 483202.03 * T) +
                        0.28 * sind(228.2 + 960400.87 * T) -
                        0.28 * sind(318.3 + 6003.18 * T) -
                        0.17 * sind(217.6 - 407332.2 * T)
    horizontal_parallax = 0.9508 + 0.0518 * cosd(134.9 + 477198.85 * T) +
                          0.0095 * cosd(259.2 - 413335.38 * T) +
                          0.0078 * cosd(235.7 + 890534.23 * T) +
                          0.0028 * cosd(269.9 + 954397.7 * T)
    ecliptic_obl = 23.439291 - 0.0130042 * T - 1.64e-7 * T^2 + 5.04e-7 * T^3
    rmoon = Rearth / sind(horizontal_parallax)
    r1 = rmoon * cosd(ecliptic_latitude) * cosd(ecliptic_longitude)
    r2 = rmoon * (cosd(ecliptic_obl) * cosd(ecliptic_latitude) *
                  sind(ecliptic_longitude) - sind(ecliptic_obl) *
                                             sind(ecliptic_latitude))
    r3 = rmoon * (sind(ecliptic_obl) * cosd(ecliptic_latitude) *
                  sind(ecliptic_longitude) + cosd(ecliptic_obl) *
                                             sind(ecliptic_latitude))
    return [r1, r2, r3]
end

function sun_pos(JD)
    JDUT1 = convert_jd(JD, :UT1)
    JDTDB = convert_jd(JD, :TDB)
    T = ((JDUT1.epoch[1] - 2451545.0) + JDUT1.epoch[2]) / 36525
    TTDB = ((JDTDB.epoch[1] - 2451545.0) + JDTDB.epoch[2]) / 36525

    ecliptic_longitude = 280.46 + 36000.771 * T
    MSUN = 357.5288 + 35999.050957 * TTDB
    ecliptic_longitude += (1.915 * sind(MSUN) + 0.02 * sind(2 * MSUN))

    rsun = 1.00014 - 0.01671 * cosd(MSUN) - 0.00014 * cosd(2 * MSUN)
    rsun *= 149597870.0
    ecliptic_obl = 23.439291 - 0.01461 * TTDB

    r1 = rsun * cosd(ecliptic_longitude)
    r2 = rsun * cosd(ecliptic_obl) * sind(ecliptic_longitude)
    r3 = rsun * sind(ecliptic_obl) * sind(ecliptic_longitude)

    return [r1, r2, r3]
end
