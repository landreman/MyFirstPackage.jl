using MyFirstPackage
using Test
using Zygote
using ForwardDiff

@testset "Optimize curve shape" begin
    nfp = 4
    curve = CurveRZFourier(
        nfp,
        [1.0, 0.05, 0],
        [0.0, 0.05, 0]
    )
    x0 = get_dofs(curve)
    nt = 15
    data = compute_curve_data(curve, nt)

    function objective(dofs)
        #set_dofs!(curve, dofs)
        new_curve = CurveRZFourier(
            nfp,
            dofs[1:3],
            [0.0; dofs[4:5]]
        )
        
        data = compute_curve_data(new_curve, nt)
        return (
            data.mean_squared_curvature 
            + 10 * (data.length - 2π) ^ 2
            + 10 * (data.integrated_torsion - 0) ^ 2
        )
        
        #return curve.rc[2]
    end

    println("Initial condition:")
    println("  length: $(data.length)")
    println("  mean squared curvature: $(data.mean_squared_curvature)")
    println("  integrated_torsion: $(data.integrated_torsion)")
    println("  objective: $(objective(x0))")

    n = 5
    ϵ = 1e-6
    finite_diff = zeros(n)
    for j in 1:n
        x = copy(x0)
        x[j] = x0[j] + ϵ
        obj_plus = objective(x)
        x[j] = x0[j] - ϵ
        obj_minus = objective(x)
        finite_diff[j] = (obj_plus - obj_minus) / (2ϵ)
    end
    println("Finite difference gradient:")
    println(finite_diff)
    println("Autodiff gradient:")
    println(gradient(objective, x0))
    #println(ForwardDiff.gradient(objective, x0))
end