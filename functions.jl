using Oscar

function witness_optimised(h, H, G, ord)
    # Assume first entry of H is a constant integer h0
    # Convert H\{h0} to Finite field mod h0

    # Divide h by h0 to get q0
    # Apply division algorithm to h w.r.t. H\{h0} to get {q1, ..., qn}

    # Lift q1, ..., qn back to original ring (ZZ[x1, ..., xn]) h = q0*h0 + q1*h1 + ... qn*hn
    # return g = q0*g0 + q1*g1 + ... qn*gn

    #If H contains multiple integers:
    #println(H[1])
    #constants = filter(Oscar.is_constant, H)
    #Hprime = filter(!Oscar.is_constant, H)

    # Assume first entry of H is a constant integer h0
    @req Oscar.is_constant(H[1]) "H[1] is not a constant polynomial"
    m = coeff(H[1], 1) # 1st coefficient of H[1] is the constant (for type purposes)

    R = parent(h)
    @req typeof(R) == Oscar.ZZMPolyRing "h is not an integer polynomial (in ZZ[x1, ..., xn])"

    S, _vars = polynomial_ring(GF(m), nvars(R)) # Create polynomial ring over finite field mod m

    phi = hom(R, S, c -> GF(m)(c), gens(S)) # homomorphism from R to S
    hbar = phi(h) # Convert h to Finite field mod h0
    Hbar = [phi(_h) for _h in H[2:end]] # Convert H\{h0} to Finite field mod h0

    q0 = div(h, H[1]) # Divide h by h0 to get q0

    u, Q, r = reduce_with_quotients_and_unit(hbar, Hbar, ordering=ord) # Apply division algorithm to h w.r.t. H\{h0} to get {q1, ..., qn} (in finite field for performane reasons)
    @req iszero(r) "Remainder is not zero"
    @req isone(u) "Unit is not one"

    phi_inv = hom(S, R, c -> lift(ZZ, c), gens(R)) # homomorphism from S to R
    Q = [phi_inv(q) for q in Q] # Lift q1, ..., qn back to original ring (ZZ[x1, ..., xn])
    Q = pushfirst!(vec(Q), q0)

    return sum([qi * gi for (qi, gi) in zip(Q, G)]) # g = q0*g0 + q1*g1 + ... qn*gn. sum(Q .* G) - (ERROR) Broadcasting for type Vector{ZZMPolyRingElem} not implemented
end

function witness(h, H, G, ord)
    u, Q, r = reduce_with_quotients_and_unit(h, H, ordering=ord)

    @req iszero(r) "Remainder is not zero"
    @req isone(u) "Unit is not one"

    return sum([qi * gi for (qi, gi) in zip(Q, G)]) # sum(Q .* G) - (ERROR) Broadcasting for type Vector{ZZMPolyRingElem} not implemented
end


function lift_custom(Hprime, ordprime, H, G, ord)
    Gprime = [witness(hprime, H, G, ord) for hprime in Hprime]
    #println(Gprimeprime)
    Gprime = reduce(Gprimeprime, ordering=ordprime) # Reduce initially?
    return Gprime
end

function flip(G, H, v, ord)
    I = ideal(H) #w inital ideal of <H>?
    ord_w = weight_ordering(w, lex(R)) # R?
    ord_wv = weight_ordering(v, ord_w)

    Hprime = standard_basis(I, ordering=ord_wv)
    Gprime = lift_custom(Hprime, ord_wv, H, G, ord)
    return (Gprime, ord_wv)
end