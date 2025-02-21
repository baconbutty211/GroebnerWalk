using Oscar

function witness(h, H, G, ord)
    # Assume first entry of H is a constant integer h0
    # Convert H\{h0} to Finite field mod h0
    # Apply division algorithm to h w.r.t. H\{h0} to get {q1, ..., qn}
    # Lift back to h_bar = q0*h0 + q1*h1 + ... qn*hn
    u, Q, r = reduce_with_quotients_and_unit(h, H, ordering=ord)
    #println(u)
    @req iszero(r) "Remainder is not zero"
    @req isone(u) "Unit is not one"
    return sum([qi * gi for (qi, gi) in zip(Q, G)])
end


function lift(Hprime, ordprime, H, G, ord)
    Gprime = [witness(hprime, H, G, ord) for hprime in Hprime]
    #println(Gprimeprime)
    Gprime = reduce(Gprimeprime, ordering=ordprime) # Reduce initially?
    return Gprime
end

function flip(G, H, v, ord)
    I = ideal(H) #w inital ideal of <H>?
    ord_w = weight_ordering(w, lex(R)) # R?
    ord_wv = weight_ordering(v, ord_w)

    Hprime = groebner_basis(I, ordering=ord_wv)
    Gprime = lift(Hprime, ord_wv, H, G, ord)
    return (Gprime, ord_wv)
end