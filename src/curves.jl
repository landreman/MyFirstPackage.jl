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

function curve_properties(curve, t0)
    #r_prime = jacobian(tt -> position_vector(curve, tt), t)[1]
    f(t) = position_vector(curve, t)
    #r_prime = ForwardDiff.derivative(f, t0)
    #r_prime_prime = ForwardDiff.derivative(t -> ForwardDiff.derivative(f, t), t0)
    #r_prime_prime_prime = ForwardDiff.derivative(t -> ForwardDiff.derivative(ForwardDiff.derivative(f, t), t), t0)
    f_prime(t) = ForwardDiff.derivative(f, t)
    f_prime_prime(t) = ForwardDiff.derivative(f_prime, t)
    f_prime_prime_prime(t) = ForwardDiff.derivative(f_prime_prime, t)

    n = length(t0)
    differential_arclength = zeros(n)
    curvature = zeros(n)
    torsion = zeros(n)
    for j in 1:n
        r_prime = f_prime(t0[j])
        r_prime_prime = f_prime_prime(t0[j])
        r_prime_prime_prime = f_prime_prime_prime(t0[j])
        
        norm_r_prime = norm(r_prime)
        differential_arclength[j] = norm_r_prime
        r_prime_cross_r_prime_prime = cross(r_prime, r_prime_prime)
        norm_r_prime_cross_r_prime_prime = norm(r_prime_cross_r_prime_prime)
        curvature[j] = (norm_r_prime_cross_r_prime_prime 
                    / (norm_r_prime * norm_r_prime * norm_r_prime))

        torsion[j] = (dot(r_prime_cross_r_prime_prime, r_prime_prime_prime) 
                    / (norm_r_prime_cross_r_prime_prime * norm_r_prime_cross_r_prime_prime))
    end

    return Dict(
            "differential_arclength" => differential_arclength,
            "curvature" => curvature,
            "torsion" => torsion,
            )
end