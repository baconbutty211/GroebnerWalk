using Oscar
R, (x, y, z, u, v) = QQ["x", "y", "z", "u", "v"]
M = gens(ideal(gens(R))^2)
I = ideal([sum(rand(Int8, length(M)) .* M) for i in 1:4])
G1 = standard_basis(I, ordering=degrevlex(R))
G2 = standard_basis(I, ordering=lex(R))
length(G1)                 # = 12  
maximum(total_degree.(G1)) # = 5
length(G2)                 # = 37
maximum(total_degree.(G2)) # = 16





function evaluate_ideal(I)
    G1 = standard_basis(I, ordering=degrevlex(R))
    println("DegRevLex - ", "length:", length(G1), ", ", "degree:", maximum(total_degree.(G1))) # = 12, 5
    G2 = standard_basis(I, ordering=lex(R))
    println("Lex - ", "length:", length(G2), ",", "degree:", maximum(total_degree.(G2))) # = 37, 16
    #    length(G1)                 # = 12
    #    maximum(total_degree.(G1)) # = 5
    #    length(G2)                 # = 37
    #    maximum(total_degree.(G2)) # = 16
end

I = ideal([x * z^2 + x * y^3 + x * y * z + y * z + x^3 + 3])
standard_basis(I, ordering=degrevlex(R))
evaluate_ideal(I)