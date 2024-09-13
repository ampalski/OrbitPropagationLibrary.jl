function mean_to_eccentric(meanAnom, e)
    eccAnomi = meanAnom
    eccAnomf = eccAnomi + (meanAnom - eccAnomi + e * sin(eccAnomi)) /
                          (1.0 - e * cos(eccAnomi))

    ctr = 1
    while abs(eccAnomf - eccAnomi) > 1e-6 && ctr < 100
        ctr += 1
        eccAnomi = eccAnomf
        eccAnomf = eccAnomi + (meanAnom - eccAnomi + e * sin(eccAnomi)) /
                              (1.0 - e * cos(eccAnomi))
    end
    return eccAnomf
end

function eccentric_to_true(eccAnom, e)
    numerator = sqrt(1 - e^2) * sin(eccAnom)
    denominator = cos(eccAnom) - e
    return atan(numerator, denominator)
end

function true_to_eccentric(trueAnom, e)
    numerator = sqrt(1 - e^2) * sin(trueAnom)
    denominator = e + cos(trueAnom)
    return atan(numerator, denominator)
end

function eccentric_to_mean(eccAnom, e)
    return eccAnom - e * sin(eccAnom)
end

function true_to_hyperbolic(trueAnom, e)
    tanhH2 = tan(trueAnom / 2) * sqrt((e - 1) / (e + 1))
    return 2 * atanh(tanhH2)
end
function hyperbolic_to_true(hypAnom, e)
    tanf2 = sqrt((e + 1) / (e - 1)) * tanh(hypAnom / 2)
    return 2 * atan(tanf2)
end
function hyperbolic_to_mean(hypAnom, e)
    return e * sinh(hypAnom) - hypAnom
end
function mean_to_hyperbolic(meanAnom, e)
    hypAnomi = meanAnom
    hypAnomf = hypAnomi + (meanAnom + hypAnomi - e * sinh(hypAnomi)) /
                          (e * cosh(hypAnomi) - 1.0)

    ctr = 1
    while abs(hypAnomf - hypAnomi) > 1e-6 && ctr < 100
        ctr += 1
        hypAnomi = hypAnomf
        hypAnomf = hypAnomi + (meanAnom + hypAnomi - e * sinh(hypAnomi)) /
                              (e * cosh(hypAnomi) - 1.0)
    end
    return hypAnomf
end

export state_to_classical
function state_to_classical(r, v)
    rn = norm(r)
    vn = norm(v)
    h = cross(r, v)
    #
    #node vector
    n = cross([0.0, 0.0, 1], h)
    e = ((vn^2 - μ / rn) * r - dot(r, v) * v) / μ
    ecc = norm(e)

    alpha = 2 / rn - vn^2 / μ

    abs(alpha) < 1e-12 && error("Not built to handle parabolic orbits")

    a = 1 / alpha
    i = acos(h[3] / norm(h))
    raan = acos(n[1] / norm(n))
    argp = acos(dot(n, e) / norm(n) / ecc)
    trueAnom = acos(clamp(dot(e, r) / ecc / rn, -1.0, 1.0))

    if n[2] < 0
        raan = 2 * pi - raan
    end
    if e[3] < 0
        argp = 2 * pi - argp
    end
    if dot(r, v) < 0
        trueAnom = 2 * pi - trueAnom
    end

    if alpha > 0
        meanAnom = eccentric_to_mean(true_to_eccentric(trueAnom, ecc), ecc)
    else
        meanAnom = hyperbolic_to_mean(true_to_hyperbolic(trueAnom, ecc), ecc)
    end
    return SA[a, ecc, i, raan, argp, meanAnom, trueAnom]
end
