#using Zygote
using ForwardDiff

mutable struct CurveRZFourier
    nfp::Integer
    rc::Vector
    zs::Vector
end

function position_vector(curve::CurveRZFourier, t)
    r = 0.0
    z = 0.0
    for j in 1:length(curve.rc)
        r += curve.rc[j] * cos(curve.nfp * (j - 1) * t)
        z += curve.zs[j] * sin(curve.nfp * (j - 1) * t)
    end
    return [r * cos(t), r * sin(t), z]
end

function get_dofs(curve::CurveRZFourier)
    return [curve.rc; curve.zs[2:end]]
end

function set_dofs(curve::CurveRZFourier, dofs)
    @assert length(dofs) % 2 == 1
    n = Integer(ceil(length(dofs) / 2))
    curve.rc = dofs[1 : n]
    curve.zs = [0.0; dofs[n + 1 : end]]
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
    t = collect(range(0, 2π / nfp, nt + 1))[1 : end - 1]
    dt = t[2] - t[1]

    #r_prime = jacobian(tt -> position_vector(curve, tt), t)[1]
    f(tt) = position_vector(curve, tt)
    #r_prime = ForwardDiff.derivative(f, t0)
    #r_prime_prime = ForwardDiff.derivative(t -> ForwardDiff.derivative(f, t), t0)
    #r_prime_prime_prime = ForwardDiff.derivative(t -> ForwardDiff.derivative(ForwardDiff.derivative(f, t), t), t0)
    f_prime(tt) = ForwardDiff.derivative(f, tt)
    f_prime_prime(tt) = ForwardDiff.derivative(f_prime, tt)
    f_prime_prime_prime(tt) = ForwardDiff.derivative(f_prime_prime, tt)

    differential_arclength = zeros(nt)
    curvature = zeros(nt)
    torsion = zeros(nt)
    for j in 1:nt
        r_prime = f_prime(t[j])
        r_prime_prime = f_prime_prime(t[j])
        r_prime_prime_prime = f_prime_prime_prime(t[j])
        
        norm_r_prime = norm(r_prime)
        differential_arclength[j] = norm_r_prime
        r_prime_cross_r_prime_prime = cross(r_prime, r_prime_prime)
        norm_r_prime_cross_r_prime_prime = norm(r_prime_cross_r_prime_prime)
        curvature[j] = (norm_r_prime_cross_r_prime_prime 
                    / (norm_r_prime * norm_r_prime * norm_r_prime))

        torsion[j] = (dot(r_prime_cross_r_prime_prime, r_prime_prime_prime) 
                    / (norm_r_prime_cross_r_prime_prime * norm_r_prime_cross_r_prime_prime))
    end

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