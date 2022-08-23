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
    end

    @testset "Test set_dofs and get_dofs" begin
        curve = CurveRZFourier(2, [2.0, 0.6, -0.3], [0.0, -1.1, 0.4])
        @test get_dofs(curve) ≈ [2.0, 0.6, -0.3, -1.1, 0.4]
        
        x = rand(3)
        set_dofs(curve, x)
        @test get_dofs(curve) == x
        
        x = rand(7)
        set_dofs(curve, x)
        @test get_dofs(curve) == x
    end
end

@testset "Test curve_properties" begin
    @testset "Check a circle" begin
        R = 2.3
        curve = CurveRZFourier(3, [R], [0])
        n = 13
        t = collect(range(0.0, 2π, length=n))
        properties = curve_properties(curve, t)
        @test properties["differential_arclength"] ≈ ones(n) * R
        @test properties["curvature"] ≈ ones(n) * 1 / R
        @test properties["torsion"] ≈ zeros(n)
    end

    @testset "Compare to QSC" begin
        # O(r^1 section 5.2)
        nfp = 2
        curve = CurveRZFourier(
            nfp,
            [1.     , 0.173  , 0.0168 , 0.00101],
            [0.      , 0.159   , 0.0165  , 0.000985]
            )
        t = collect(range(0, π, length=16))[1:end - 1]
        properties = curve_properties(curve, t)

        @test properties["differential_arclength"] ≈ [1.25301966, 1.23279065, 1.17766599, 1.10202631, 1.02332187,
            0.95601355, 0.90857199, 0.88440396, 0.88440396, 0.90857199,
            0.95601355, 1.02332187, 1.10202631, 1.17766599, 1.23279065]

        @test properties["curvature"] ≈ [1.39355975, 1.36904437, 1.29001429, 1.14849644, 0.95518636,
            0.75252884, 0.59242666, 0.5048237 , 0.5048237 , 0.59242666,
            0.75252884, 0.95518636, 1.14849644, 1.29001429, 1.36904437]

        @test properties["torsion"] ≈ [-0.40604021, -0.42808217, -0.46686814, -0.4364834 , -0.20478637,
            0.31002864,  0.99703528,  1.58019348,  1.58019348,  0.99703528,
            0.31002864, -0.20478637, -0.4364834 , -0.46686814, -0.42808217]
    end
end