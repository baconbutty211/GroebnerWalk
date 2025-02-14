using Oscar

function witness(h, H, G, ord)
    u, Q, r = reduce_with_quotients_and_unit(h, H, ordering=ord)
    #println(u)
    @req iszero(r) "Remainder is not zero"
    @req isone(u) "Unit is not one"
    return sum([qi * gi for (qi, gi) in zip(Q, G)])
end

function lift(H', ord', H, G, ord)
    G'' = [witness(h', H, G, ord) for h' in H']
    #println(G'')
    G' = reduce(G'', ordering=ord') # Reduce initially?
    return G'
end

function flip(G, H, v, ord)
    I = ideal(H) #w inital ideal of <H>?
    ord_w = weight_ordering(w, lex(R)) # R?
    ord_wv = weight_ordering(v, ord_w)

    H' = groebner_basis(I, ordering=ord_wv)
    G' = lift(H', ord_wv, H, G, ord)
end