include("./GroebnerWalk.jl")
include("./xHomoGWalk.jl")
include("./functions.jl")
include("./initialReduction.jl")

function witness_optimised_test()
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

    print(G[1] == witness_optimised(h_1, H, G, ord))
    print(G[2] == witness_optimised(h_2, H, G, ord))
end

function witness_optimised_test2()
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
end

function is_reduced_test()
    R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
    G = [x^2 - y, y^2 + 1]

    println(isReduced(G, lex(R))) # Check if G is reduced w.r.t. lex ordering
end

function is_initially_reduced_test()
    # Test isInitiallyReduced function
    R, (t, x1, x2, x3) = polynomial_ring(ZZ, ["t", "x1", "x2", "x3"])

    ord = weight_ordering([-1, 1, 1, 1], lex(R))

    g1 = (1 - 3t^2 + t^4 - t^5 + t^6 + t^7) * x1^2
    g2 = (1 - t^2) * x2^2 + (t + t^2 + t^3) * x3^2
    g3 = (t^3 - 2t^5 - t^7 - t^8) * x3^2
    G = [g1, g2, g3]

    println(isInitiallyReduced(G, ord)) # Check if G is initially reduced w.r.t <G, 2-t> under lex ordering
end

function same_degree_initially_reduced_test()
    R, (t, x1, x2, x3) = polynomial_ring(ZZ, ["t", "x1", "x2", "x3"])

    g1 = x1^2 + t * x2^2 - t^2 * x3^2
    g2 = t * x1^2 + x2^2 + (t * x3^2 + t^2 * x3^2)
    g3 = t^4 * x1^2 + (t^4 * x2^2 + t^5 * x2^2) + t^3 * x3^2
    G = [g1, g2, g3]

    ord = weight_ordering([-1, 1, 1, 1], lex(R))

    G_init_red_actual = sameDegreeReduce(G, 2, ord) # Apply alg 4.2 to initially reduce G w.r.t <G, 2-t> under ord

    g1_expected = (1 - 5t^2 + 3t^4 - t^5 + t^6 + t^7) * x1^2
    g2_expected = (1 - t^2) * x2^2 + (t + t^2 + t^3) * x3^2
    g3_expected = (t^3 - 2t^5 - t^7 - t^8) * x3^2
    G_init_red_expected = [g1_expected, g2_expected, g3_expected]

    println(isInitiallyReduced(G_init_red_actual, ord)) # Check if the initially reduced G is equal to the expected result
    println(G_init_red_actual == G_init_red_expected) # Check if the initially reduced G is equal to the expected result
end

function all_at_once_initially_reduced_test()
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
end

function initially_reduce_test()
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
end

function flip_test()
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
end


# (TODO) Test minimise()
# (TODO) Test isMinimised()