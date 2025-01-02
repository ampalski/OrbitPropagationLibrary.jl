# Will need to pass options, if not the whole input, through p
function force_model(x, p, t)
    accel = zeros(3)
    jd0 = p[1]
    opts = p[2]
    jd = JDate(SA[jd0.epoch[1], jd0.epoch[2]+t/86400.0], jd0.system)

    r = norm(x[1:3])
    a_twobody = -μ / r^3 * x[1:3]
    accel += a_twobody

    if opts.third_body_moon

    end

    return SA[x[4:6]...; accel...]
end
