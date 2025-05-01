using Oscar
include("./xHomoGWalk.jl")

R, (t, x, y) = polynomial_ring(ZZ, ["t", "x", "y"])
u = [-1, 1, 1] # Initial Weight vector for the ordering
ord_t_lex = matrix_ordering(R, [0 1 0; 0 0 1; -1 0 0]) # t-local lexicographic ordering
ord_u = weight_ordering(u, ord_t_lex) # Create a new ordering with u as the weight vector

I = ideal([2 - t, x * y^2 - t^2 * y^3, x^2 - t^3 * y^2])
G = [2 - t, x * y^2 - t^2 * y^3, x^2 - t^3 * y^2, t^3 * y^4] # Grobner basis of I w.r.t. >_u


w = [-4, 1, 7] # Weight vector on the boundary of Grobner cone for the ordering
ord = weight_ordering(w, ord_t_lex)
v = [3, 5, 1] # Outward pointing normal vector in positive orthant

H = initial(collect(G), ord, ZZ.(w))

G_expected = [2 - t, x * y^2 - t^2 * y^3, t^3 * y^2 - x^2, x^3 - t^5 * y^3] # Expected result of the flip

(G_prime, ord_prime) = flip(G, H, v, ord)
@req isInitiallyReduced(G_prime, ord) "G' is not initially reduced w.r.t. ord"
@req G_pime[1] == G_expected[1] "G'[1]=$(G_pime[1]) is not equal to the expected result $(G_expected[1])"
@req G_pime[2] == G_expected[2] "G'[2]=$(G_pime[2]) is not equal to the expected result $(G_expected[2])"
@req G_pime[3] == G_expected[3] "G'[3]=$(G_pime[3]) is not equal to the expected result $(G_expected[3])"
@req G_pime[4] == G_expected[4] "G'[4]=$(G_pime[4]) is not equal to the expected result $(G_expected[4])"

(G_prime, ord_prime) = flip_optimised(G, H, v, ord)
@req G_pime[1] == G_expected[1] "G'[1]=$(G_pime[1]) is not equal to the expected result $(G_expected[1])"
@req G_pime[2] == G_expected[2] "G'[2]=$(G_pime[2]) is not equal to the expected result $(G_expected[2])"
@req G_pime[3] == G_expected[3] "G'[3]=$(G_pime[3]) is not equal to the expected result $(G_expected[3])"
@req G_pime[4] == G_expected[4] "G'[4]=$(G_pime[4]) is not equal to the expected result $(G_expected[4])"