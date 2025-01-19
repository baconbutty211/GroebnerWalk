include("GroebnerWalk.jl")

R, (t, x, y, z) = polynomial_ring(ZZ, ["t", "x", "y", "z"])
w = [-1, 3, 3, 3]
w_x = [0, 1, 1, 1]
v = -w + 4 * w_x
ord = weight_ordering(v, lex(R))

G = [x - t^3 * x + t^3 * z - t^4 * z, y - t^3 * y + t^2 * z - t^4 * z]
H = [x, y + t^2 * z]
h_1 = x
h_2 = y + t^2 * z

u, Q, r = reduce_with_quotients_and_unit(h_1, H, ordering=ord)
println(u)

print(G[1] == witness(h_1, H, G, ord))
print(G[2] == witness(h_2, H, G, ord))