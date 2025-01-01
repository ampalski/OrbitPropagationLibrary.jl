# Will need to pass options, if not the whole input, through p
function force_model(x, p, t)
    accel = zeros(3)

    r = norm(x[1:3])
    a_twobody = -μ / r^3 * x[1:3]
    accel += a_twobody

    return SA[x[4:6]...; accel...]
end
