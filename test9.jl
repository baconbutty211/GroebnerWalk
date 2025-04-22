include("./xHomoGWalk.jl")
include("./functions.jl")

R, (t, x1, x2, x3, x4) = polynomial_ring(ZZ, ["t", "x1", "x2", "x3", "x4"])
w = [-1, 1, 11, 3, 19]
ord = weight_ordering(w, lex(R))

w_x = [0, 1, 1, 1, 1]
v = -w + 20 * w_x # Change v to be an outward pointing normal vector

I = ideal([3 − t, 2 * x1^2 + 3 * x1 * x2 + 24 * x3 * x4, 8 * x1^3 + x2 * x3 * x4 + 18 * x3^2 * x4])
G = standard_basis(I, ordering=ord)
H = initial(collect(G), ord, ZZ.(w))

G_prime = flip(G, H, v, ord)
println(isInitiallyReduced(G_prime, ord))