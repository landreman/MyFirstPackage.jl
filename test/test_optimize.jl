using MyFirstPackage
using Test
using Zygote
using ForwardDiff

@testset "Optimize curve shape" begin
    curve = CurveRZFourier(
        4,
        [1.0, 0.05, 0],
        [0.0, 0.05, 0]
    )
    nt = 15
    data = compute_curve_data(curve, nt)

    function objective(dofs)
        set_dofs(curve, dofs)
        
        data = compute_curve_data(curve, nt)
        return (
            data.mean_squared_curvature 
            + 10 * (data.length - 2π) ^ 2
            + 10 * (data.integrated_torsion - 0) ^ 2
        )
        
        #return curve.rc[2]
    end

    x0 = get_dofs(curve)
    println("Initial condition:")
    println("  length: $(data.length)")
    println("  mean squared curvature: $(data.mean_squared_curvature)")
    println("  integrated_torsion: $(data.integrated_torsion)")
    println("  objective: $(objective(x0))")

    #println(gradient(objective, x0))
    println(ForwardDiff.gradient(objective, x0))
end