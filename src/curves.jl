using Zygote

struct CurveRZFourier
    #mutable struct CurveRZFourier
    nfp::Int
    rc::Vector{Float64}
    zs::Vector{Float64}
end

"""
    position_vector(curve::CurveRZFourier, t)

Compute the position vector along with its first 3 derivatives with respect to
the curve parameter t.

We could use automatic differentiation to get the derivatives. But this would
introduce a second level of automatic differentiation besides getting the
overall derivative of the objective function, which doesn't seem to work. So
here we just evaluate the derivatives analytically, which is not bad.
"""
function position_vector(curve::CurveRZFourier, t)
    r = [0.0, 0.0, 0.0]
    drdt = [0.0, 0.0, 0.0]
    d2rdt2 = [0.0, 0.0, 0.0]
    d3rdt3 = [0.0, 0.0, 0.0]
    cost = cos(t)
    sint = sin(t)
    for j in 1:length(curve.rc)
        n = curve.nfp * (j - 1)
        cosnt = cos(n * t)
        sinnt = sin(n * t)

        r += [
            curve.rc[j] * cosnt * cost,
            curve.rc[j] * cosnt * sint,
            curve.zs[j] * sinnt
        ]

        drdt += [
            curve.rc[j] * (-cosnt * sint - n * sinnt * cost),
            curve.rc[j] * (cosnt * cost - n * sinnt * sint),
            curve.zs[j] * n * cosnt
        ]

        d2rdt2 += [
            curve.rc[j] * (-cosnt * cost + 2 * n * sinnt * sint - n * n * cosnt * cost),
            curve.rc[j] * (-cosnt * sint - 2 * n * sinnt * cost - n * n * cosnt * sint),
            curve.zs[j] * (-n * n * sinnt)
        ]

        d3rdt3 += [
            curve.rc[j] * (cosnt * sint + 3 * n * sinnt * cost + 3 * n * n * cosnt * sint + n * n * n * sinnt * cost),
            curve.rc[j] * (-cosnt * cost + 3 * n * sinnt * sint - 3 * n * n * cosnt * cost + n * n * n * sinnt * sint),
            curve.zs[j] * (-n * n * n * cosnt)
        ]
    end
    return r, drdt, d2rdt2, d3rdt3
end

function get_dofs(curve::CurveRZFourier)
    return [curve.rc; curve.zs[2:end]]
end

function set_dofs!(curve::CurveRZFourier, dofs)
    @assert length(dofs) % 2 == 1
    n = Integer(ceil(length(dofs) / 2))
    curve.rc[:] = dofs[1 : n]
    curve.zs[:] = [0.0; dofs[n + 1 : end]]
end

dot(v1, v2) = v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]

norm(v) = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3])

function cross(v1, v2)
    return [v1[2] * v2[3] - v1[3] * v2[2],
            v1[3] * v2[1] - v1[1] * v2[3],
            v1[1] * v2[2] - v1[2] * v2[1]]
end

mutable struct curve_data
    nfp
    nt
    t
    dt
    differential_arclength
    curvature
    torsion
    length
    integrated_torsion
    mean_squared_curvature
end

curve_data() = curve_data(
    1, 1, [], 0,  # nfp, nt, t, dt
    [], [], [],  # differential_arclength, curvature, torsion
    0, 0, 0  # length, integrated_torsion, mean_squared_curvature
    )

function compute_curve_data(curve, nt)
    nfp = curve.nfp
    #t = collect(range(0, 2π / nfp, length = nt + 1))[1 : end - 1]  # Zygote isn't able to differentiate through this for some reason.
    t = [j * 2π / (nt * nfp) for j in 0:(nt - 1)]
    dt = t[2] - t[1]

    "=
    f(tt) = position_vector(curve, tt)
    f_prime(tt) = ForwardDiff.derivative(f, tt)
    f_prime_prime(tt) = ForwardDiff.derivative(f_prime, tt)
    f_prime_prime_prime(tt) = ForwardDiff.derivative(f_prime_prime, tt)
    ="

    #differential_arclength = zeros(nt)
    #curvature = zeros(nt)
    #torsion = zeros(nt)
    differential_arclength = Zygote.Buffer(t)
    curvature = Zygote.Buffer(t)
    torsion = Zygote.Buffer(t)
    for j in 1:nt
        "=
        r_prime = f_prime(t[j])
        r_prime_prime = f_prime_prime(t[j])
        r_prime_prime_prime = f_prime_prime_prime(t[j])
        ="
        r, r_prime, r_prime_prime, r_prime_prime_prime = position_vector(curve, t[j])

        norm_r_prime = norm(r_prime)
        differential_arclength[j] = norm_r_prime
        r_prime_cross_r_prime_prime = cross(r_prime, r_prime_prime)
        norm_r_prime_cross_r_prime_prime = norm(r_prime_cross_r_prime_prime)
        curvature[j] = (norm_r_prime_cross_r_prime_prime 
                    / (norm_r_prime * norm_r_prime * norm_r_prime))

        torsion[j] = (dot(r_prime_cross_r_prime_prime, r_prime_prime_prime) 
                    / (norm_r_prime_cross_r_prime_prime * norm_r_prime_cross_r_prime_prime))
    end

    differential_arclength = copy(differential_arclength)
    curvature = copy(curvature)
    torsion = copy(torsion)

    length = sum(differential_arclength) * dt * nfp
    integrated_torsion = sum(differential_arclength .* torsion) * dt * nfp
    mean_squared_curvature = (sum(differential_arclength .* curvature .* curvature)
        * dt * nfp / length)

    # Return the results in a struct
    data = curve_data()

    data.nfp = curve.nfp
    data.nt = nt
    data.t = t
    data.dt = dt

    data.differential_arclength = differential_arclength
    data.curvature = curvature
    data.torsion = torsion

    data.length = length
    data.integrated_torsion = integrated_torsion
    data.mean_squared_curvature = mean_squared_curvature
    
    return data
end