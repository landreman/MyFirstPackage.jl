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
    return r * cos(t), r * sin(t), z
end

function get_dofs(curve::CurveRZFourier)
    return [curve.rc; curve.zs[2:end]]
end