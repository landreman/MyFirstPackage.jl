module MyFirstPackage

#export fun1, fun2, spectral_diff_matrix
export fun1, fun2
export CurveRZFourier, position_vector, get_dofs, set_dofs, curve_properties

include("functions.jl")
include("spectral_diff_matrix.jl")
include("curves.jl")

function fun1(x, y)
    return 2 * x + y
end

end
