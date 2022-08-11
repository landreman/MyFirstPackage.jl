module MyFirstPackage

export fun1, fun2, spectral_diff_matrix

include("functions.jl")
include("spectral_diff_matrix.jl")

function fun1(x, y)
    return 2 * x + y
end

end
