using Oscar

function witness(h, H, G, ord)
    u, Q, r = reduce_with_quotients_and_unit(h, H, ordering=ord)
    #println(u)
    @req iszero(r) "Remainder is not zero"
    @req isone(u) "Unit is not one"
    return sum([qi * gi for (qi, gi) in zip(Q, G)])
end

function lift()

end

function flip()

end