using MyFirstPackage
using Test

@testset "Test CurveRZFourier" begin
    @testset "Check position vector for a specific case" begin
        nfp = 5
        curve = CurveRZFourier(nfp, [2.0, 0.6, -0.3], [0.0, -1.1, 0.4])
        for t in range(0.0, 2π, length=13)
            x1, y1, z1 = position_vector(curve, t)
            r = 2.0 + 0.6 * cos(nfp * t) - 0.3 * cos(2 * nfp * t)
            z2 = -1.1 * sin(nfp * t) + 0.4 * sin(2 * nfp * t)
            x2 = r * cos(t)
            y2 = r * sin(t)
            @test x1 ≈ x2
            @test y1 ≈ y2
            @test z1 ≈ z2
        end
        @test get_dofs(curve) ≈ [2.0, 0.6, -0.3, -1.1, 0.4]
    end
end