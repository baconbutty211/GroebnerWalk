# Test isInitiallyReduced function
include("./initialReduction.jl")
R, (t, x1, x2, x3) = polynomial_ring(ZZ, ["t", "x1", "x2", "x3"])

ord = weight_ordering([-1, 1, 1, 1], lex(R))

g1 = (1 - 3t^2 + t^4 - t^5 + t^6 + t^7) * x1^2
g2 = (1 - t^2) * x2^2 + (t + t^2 + t^3) * x3^2
g3 = (t^3 - 2t^5 - t^7 - t^8) * x3^2
G = [g1, g2, g3]

println(isInitiallyReduced(G, ord)) # Check if G is initially reduced w.r.t <G, 2-t> under lex ordering