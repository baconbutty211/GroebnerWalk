using Oscar


function lift_to_coefficient_ring(R, f)
    C = coefficient_ring(R)
    if iszero(f)
        return zero(C)
    end

    coeffs = [lift(C, c) for c in coefficients(f)]
    monos = [prod(gens(R) .^ degrees(m)) for m in monomials(f)]
    fbar = sum(coeffs .* monos)

    #println("coefficients type: ", typeof(coeffs))
    #println("monomials type: ", typeof(monos))
    #println("lifted polynomial: ", typeof(fbar))

    return fbar
end

function witness_optimised(h, H, G, ord)
    # Assume first entry of H is a constant integer h0
    # Convert H\{h0} to Finite field mod h0
    # R,(x,y) = polynomial_ring(ZZ,["x","y"])
    # S,(x,y) = polynomial_ring(GF(3),["x","y"])
    # phi = hom(R,S,c->GF(3)(c),[x,y]) # homomorphism from R to S
    # Apply division algorithm to h w.r.t. H\{h0} to get {q1, ..., qn}
    # Lift back to h_bar = q0*h0 + q1*h1 + ... qn*hn

    #println(H[1])
    #constants = filter(Oscar.is_constant, H)
    #Hprime = filter(!Oscar.is_constant, H)

    @req Oscar.is_constant(H[1]) "H[1] is not a constant polynomial"
    R = parent(h)

    c = coeff(H[1], 1) # 1st coefficient of H[1] is the constant (for type purposes)
    hbar = change_coefficient_ring(GF(c), h)
    Hbar = [change_coefficient_ring(GF(c), hprime) for hprime in H[2:end]] # Convert to finite field mod c

    q0 = div(h, H[1])
    u, Q, r = reduce_with_quotients_and_unit(hbar, Hbar, ordering=ord) # Division algorithm (in finite field for performane reasons)

    Q = [lift_to_coefficient_ring(R, q) for q in Q] # Conversion back to original ring (ZZ[x1, ..., xn])
    Q = pushfirst!(vec(Q), q0)

    #println(typeof(r), r)
    #println(typeof(u), u)
    @req iszero(r) "Remainder is not zero"
    @req isone(u) "Unit is not one"

    return sum([qi * gi for (qi, gi) in zip(Q, G)]) # sum(Q .* G) - (ERROR) Broadcasting for type Vector{ZZMPolyRingElem} not implemented
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