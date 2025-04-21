include("./initialReduction.jl")

R, (t, x1, x2, x3) = polynomial_ring(ZZ, ["t", "x1", "x2", "x3"])
ord = weight_ordering([-1, 1, 1, 1], lex(R))

h1 = x1^2 + t * x2^2 - t^2 * x3^2
h2 = t * x1^2 + x2^2 + (t * x3^2 + t^2 * x3^2)
h3 = t^4 * x1^2 + (t^4 * x2^2 + t^5 * x2^2) + t^3 * x3^2
H = [h1, h2, h3]

g1 = x1 + t * x2 - t^2 * x3
g2 = t * x1 + x2 + (t + t^2) * x3
g3 = t^4 * x1 + (t^4 + t^5) * x2 + t^3 * x3
G = [g1, g2, g3]

H_init_red_actual = allAtOnceReduce(G, H, 2, ord) # Apply alg 4.2 to initially reduce G w.r.t <G, 2-t> under ord
println(isInitiallyReduced(H_init_red_actual, ord)) # Check if the initially reduced G is equal to the expected result


#g1_expected = (1 - 5t^2 + 3t^4 - t^5 + t^6 + t^7) * x1^2
#g2_expected = (1 - t^2) * x2^2 + (t + t^2 + t^3) * x3^2
#g3_expected = (t^3 - 2t^5 - t^7 - t^8) * x3^2
#G_init_red_expected = [g1_expected, g2_expected, g3_expected]
#println(G_init_red_actual == G_init_red_expected) # Check if the initially reduced G is equal to the expected result