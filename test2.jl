include("./xHomoGWalk.jl")

R, (t, x1, x2, x3, x4) = polynomial_ring(ZZ, ["t", "x1", "x2", "x3", "x4"])
w = [-1, 1, 11, 3, 19]
w_x = [0, 1, 1, 1, 1]
v = -w + 20 * w_x
#ord = weight_ordering(v, lex(R))
ord = weight_ordering(w, lex(R))

I = ideal([3 − t, 2 * x1^2 + 3 * x1 * x2 + 24 * x3 * x4, 8 * x1^3 + x2 * x3 * x4 + 18 * x3^2 * x4])
G = standard_basis(I, ordering=ord)
H = initial(collect(G), ord, ZZ.(w))

println(G[1] == witness(H[1], H, G, ord))
println(G[2] == witness(H[2], H, G, ord))
println(G[3] == witness(H[3], H, G, ord))
println(G[4] == witness(H[4], H, G, ord))

println(G[1] == witness_optimised(H[1], H, G, ord))
println(G[2] == witness_optimised(H[2], H, G, ord))
println(G[3] == witness_optimised(H[3], H, G, ord))
println(G[4] == witness_optimised(H[4], H, G, ord))