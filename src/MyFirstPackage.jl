module MyFirstPackage

export fun1, fun2

include("functions.jl")

function fun1(x, y)
    return 2 * x + y
end

end
