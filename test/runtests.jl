using MyFirstPackage
using Test

@testset "Test spectral_diff_matrix" begin
 
    @testset "Test matrix for n = 2" begin
        @test MyFirstPackage.spectral_diff_matrix(2) ≈ zeros(2, 2) atol=1e-13
    end

    @testset "Test matrix for n = 3" begin
        x = 0.577350269189626
        D = [ 0;  x;  -x;;
             -x;  0;   x;;
              x; -x;   0]
        @test MyFirstPackage.spectral_diff_matrix(3) ≈ D atol=1e-13
    end

    @testset "Test matrix for n = 4" begin
        D = [0; 0.5; 0; -0.5;;
            -0.5; 0; 0.5; 0;;
            0; -0.5; 0; 0.5;;
            0.5; 0;-0.5; 0]
        @test MyFirstPackage.spectral_diff_matrix(4) ≈ D atol=1e-13
    end

    @testset "Test matrix for n = 5" begin
        e = 0.85065080835204
        f = 0.525731112119134
        D = [0; e; -f; f; -e;;
            -e; 0; e; -f; f;;
            f; -e; 0; e; -f;;
            -f; f; -e; 0; e;;
            e; -f; f; -e; 0]
        @test MyFirstPackage.spectral_diff_matrix(5) ≈ D atol=1e-13
    end

end
