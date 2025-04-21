# Test isReduced function
include("functions.jl")

R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
G = [x^2 - y, y^2 + 1]

println(isReduced(G, lex(R))) # Check if G is reduced w.r.t. lex ordering