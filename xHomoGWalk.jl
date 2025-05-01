using Oscar
include("./initialReduction.jl")

function witness_optimised(h, H, G, ord)
    # Assume first entry of H is a constant integer h0
    # Convert H\{h0} to Finite field mod h0

    # Divide h by h0 to get q0
    # Apply division algorithm to h w.r.t. H\{h0} to get {q1, ..., qn}

    # Lift q1, ..., qn back to original ring (ZZ[x1, ..., xn]) h = q0*h0 + q1*h1 + ... qn*hn
    # return g = q0*g0 + q1*g1 + ... qn*gn

    # Assume first entry of H is a constant integer h0
    @req Oscar.is_constant(H[1]) "H[1] is not a constant polynomial"
    p = coeff(H[1], 1) # 1st coefficient of H[1] is the constant (for type purposes)

    R = parent(h)
    @req typeof(R) == Oscar.ZZMPolyRing "h is not an integer polynomial (in ZZ[x1, ..., xn])"

    S, _vars = polynomial_ring(GF(p), nvars(R)) # Create polynomial ring over finite field mod h0

    phi = hom(R, S, c -> GF(p)(c), gens(S)) # homomorphism from R to S
    hbar = phi(h) # Convert h to Finite field mod h0
    Hbar = [phi(_h) for _h in H[2:end]] # Convert H\{h0} to Finite field mod h0

    q0 = div(h, H[1]) # Divide h by h0 to get q0

    u, Q, r = reduce_with_quotients_and_unit(hbar, Hbar, ordering=ord) # Apply division algorithm to h w.r.t. H\{h0} to get {q1, ..., qn} (in finite field for performance reasons)
    @req iszero(r) "Remainder is not zero"
    @req isone(u) "Unit is not one"

    phi_inv = hom(S, R, c -> lift(ZZ, c), gens(R)) # homomorphism from S to R
    Q = [phi_inv(q) for q in Q] # Lift q1, ..., qn back to original ring (ZZ[x1, ..., xn])
    Q = pushfirst!(vec(Q), q0)

    f = sum([qi * gi for (qi, gi) in zip(Q, G)]) # g = q0*g0 + q1*g1 + ... qn*gn. sum(Q .* G) - (ERROR) Broadcasting for type Vector{ZZMPolyRingElem} not implemented
    @req h == initial(f, ord, matrix(ord)[1, :]) "Initial form of Witness is not equal to h"
    return f
end

function witness(h, H, G, ord)
    u, Q, r = reduce_with_quotients_and_unit(h, H, ordering=ord)

    @req iszero(r) "Remainder is not zero"
    @req isone(u) "Unit is not one"

    f = sum([qi * gi for (qi, gi) in zip(Q, G)]) # sum(Q .* G) - (ERROR) Broadcasting for type Vector{ZZMPolyRingElem} not implemented
    @req h == initial(f, ord, matrix(ord)[1, :]) "Initial form of Witness is not equal to h"
    return f
end


function lift_optimised(Hprime, ordprime, H, ord, G)
    Gprimeprime = [witness_optimised(hprime, H, G, ord) for hprime in Hprime]
    #println(Gprimeprime)
    # Assume h_0 = p 
    p = coeff(H[1], 1)
    Gprime = initially_reduce(Gprimeprime, p, ordprime) # No method for initial reduction is implemented in Oscar, See Algorithm 4.7 in https://arxiv.org/abs/1512.02662
    Gprime = Gprimeprime
    return Gprime
end
function lift_custom(Hprime, ordprime, H, ord, G)
    Gprimeprime = [witness(hprime, H, G, ord) for hprime in Hprime]
    #println(Gprimeprime)
    # Assume h_0 = p 
    p = coeff(H[1], 1)
    Gprime = initially_reduce(Gprimeprime, p, ordprime) # No method for initial reduction is implemented in Oscar, See Algorithm 4.7 in https://arxiv.org/abs/1512.02662
    Gprime = Gprimeprime
    return Gprime
end

function flip_optimised(G, H, v, ord_w)
    @req typeof(ord_w) <: MonomialOrdering "ord_w is not a weight ordering"

    w = matrix(ord_w)[1, :] # w is the weight vector of the ordering
    @req length(w) == nvars(parent(H[1])) "Length of weight vector w is not equal to number of variables in the polynomial ring"
    @req length(w) == length(v) "Length of weight vector w must be equal to v"
    @req w[1] < 0 "Weight vector w must be negative in the first entry"

    @req length(G) == length(H) "Length of G and H must be equal"
    @req H == initial(collect(G), ord_w, ZZ.(w)) "H is not the initial form of G w.r.t. w"

    I = ideal(H)

    R = parent(G[1])
    ord_v = weight_ordering(v, lex(R)) # Create a new ordering with v as the weight vector
    ord_wv = weight_ordering(w, ord_v) # Create a new ordering with w as the weight vector and ord_v as the tie-breaker

    Hprime = standard_basis(I, ordering=ord_wv)
    Gprime = lift_optimised(Hprime, ord_wv, H, ord, G)
    return (Gprime, ord_wv)
end
function flip(G, H, v, ord_w)
    @req typeof(ord_w) <: MonomialOrdering "ord_w is not a weight ordering"

    w = Int.(matrix(ord_w)[1, :]) # w is the weight vector of the ordering

    @req typeof(w) == Vector{Int} "w is not a vector of integers"
    @req length(w) == nvars(parent(H[1])) "Length of weight vector w is not equal to number of variables in the polynomial ring"
    @req length(w) == length(v) "Length of weight vector w must be equal to v"
    @req w[1] < 0 "Weight vector w must be negative in the first entry"

    @req length(G) == length(H) "Length of G and H must be equal"
    @req H == initial(collect(G), ord_w, ZZ.(w)) "H is not the initial form of G w.r.t. w"

    I = ideal(H)

    R = parent(G[1])
    ord_v = weight_ordering(v, lex(R)) # Create a new ordering with v as the weight vector
    ord_wv = weight_ordering(w, ord_v) # Create a new ordering with w as the weight vector and ord_v as the tie-breaker

    Hprime = standard_basis(I, ordering=ord_wv)
    Gprime = lift_custom(Hprime, ord_wv, H, ord, G)
    return (Gprime, ord_wv)
end