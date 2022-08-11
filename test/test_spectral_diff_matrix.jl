using MyFirstPackage
using Test

@testset "Test spectral_diff_matrix" begin
 
    @testset "Test matrix for n = 2" begin
        @test MyFirstPackage.spectral_diff_matrix(2) ≈ zeros(2, 2) atol=1e-13
    end

    @testset "Test matrix for n = 3" begin
        x = 0.577350269189626
        D = [ 0  x -x
             -x  0  x
              x -x  0]
        @test MyFirstPackage.spectral_diff_matrix(3) ≈ D atol=1e-13
    end

    @testset "Test matrix for n = 4" begin
        D = [   0  0.5    0 -0.5
             -0.5    0  0.5    0
                0 -0.5    0  0.5
              0.5    0 -0.5    0]
        @test MyFirstPackage.spectral_diff_matrix(4) ≈ D atol=1e-13
    end

    @testset "Test matrix for n = 5" begin
        e = 0.85065080835204
        f = 0.525731112119134
        D = [ 0  e -f  f -e
             -e  0  e -f  f
              f -e  0  e -f
             -f  f -e  0  e
              e -f  f -e  0]
        @test MyFirstPackage.spectral_diff_matrix(5) ≈ D atol=1e-13
    end

    @testset "Test matrix for n = 2 with shifted xmin and xmax" begin
        @test MyFirstPackage.spectral_diff_matrix(2, xmin=-2.1, xmax=3.7) ≈ zeros(2, 2) atol=1e-13
    end

    @testset "Test matrix for n = 3 with shifted xmin and xmax" begin
        x = 0.625448056632489
        D = [ 0  x -x
             -x  0  x
              x -x  0]
        @test MyFirstPackage.spectral_diff_matrix(3, xmin=-2.1, xmax=3.7) ≈ D atol=1e-13
    end

    @testset "Test matrix for n = 4 with shifted xmin and xmax" begin
        x = 0.541653905791344
        D = [ 0  x  0 -x
             -x  0  x  0
              0 -x  0  x
              x  0 -x  0]
        @test MyFirstPackage.spectral_diff_matrix(4, xmin=-2.1, xmax=3.7) ≈ D atol=1e-13
    end

    @testset "Test matrix for n = 5 with shifted xmin and xmax" begin
        e = 0.921516665616892
        f = 0.569528620550711
        D = [ 0  e -f  f -e
             -e  0  e -f  f
              f -e  0  e -f
             -f  f -e  0  e
              e -f  f -e  0]
        @test MyFirstPackage.spectral_diff_matrix(5, xmin=-2.1, xmax=3.7) ≈ D atol=1e-13
    end

    @testset "Confirm that the derivative of a sin/cos is computed exactly." begin
        xmin = -3.4
        xmax = -0.7
        L = xmax - xmin
        for nphi in 11:21
            D = spectral_diff_matrix(nphi, xmin=xmin, xmax=xmax)
            for n in range(0, Int(floor(nphi / 2)) - 1)
                for phase in [0, 0.3]
                    phi = range(xmin, xmax, nphi + 1)[1:nphi]
                    x = @. sin(n * phi * 2 * π / L + phase)
                    dx = @. (n * 2 * π / L) * cos(n * phi * 2 * π / L + phase)
                    @test D * x ≈ dx rtol=1e-12  atol = 1e-12
                end
            end
        end
    end
end
