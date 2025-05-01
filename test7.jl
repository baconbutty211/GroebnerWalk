include("./initialReduction.jl")
include("./functions.jl")

R, (t, x1, x2, x3) = polynomial_ring(ZZ, ["t", "x1", "x2", "x3"])
ord = weight_ordering([-1, 1, 1, 1], lex(R))

I = ideal([3 - t,
    2 * x1^2 + 3 * x1 * x2 + 24 * x3^2,
    8 * x1^3 + x2 * x3^2 + 18 * x3^3])
F = gens(I)

G_init_red_actual = initialReduce(F, 3, ord) # Apply alg 4.7 to get initially reduced standard basis of I w.r.t. ord

println(isInitiallyReduced(G_init_red_actual, ord)) # Check if the initially reduced G is equal to the expected result
# Check if G_init_red_actual is standard basis of I
# Check if G_init_red_actual is expected result