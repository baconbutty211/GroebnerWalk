using Oscar
include("./functions.jl")

function ptReduce(g, p, ord)
    # Algorithm 4.1 from https://arxiv.org/abs/1512.02662
    # Input: g ∈ R x-homogeneous, ord a t-local monomial ordering on R
    # Output: g' ∈ R x-homogeneous and initially reduced wrt p-t under > (no term of tail(g') is divisible by p) with <p-t, g'> = <p-t, g> and LT(g')= LT(g)

    g_gamma = collectLeadingXTerms(g, ord) # Collect the leading x-terms of g w.r.t. ord

    g_prime = g_gamma # Initialize g' to g_gamma
    g_primeprime = g - g_gamma # Initialize g'' to g - g_gamma

    while !iszero(g_primeprime)
        g_primeprime_gamma = lt(g_primeprime, ord) # Get the leading term of g'' w.r.t. ord
        lead_coeff = coeff(g_primeprime_gamma, 1)

        l = 0 # Initialize max power of p in g_primeprime_gamma to 0
        while div(lead_coeff, p) > 0 # Check if lead_coeff is divisible by p
            lead_coeff = div(lead_coeff, p) # Divide lead_coeff by p
            l += 1 # Increment max power of p in g_primeprime_gamma
            #println("l = ", l, " lead_coeff = ", lead_coeff)
        end

        g_primeprime = g_primeprime - (mod(coeff(g_primeprime_gamma, 1), p) * (p^l - t^l))
        g_primeprime_gamma = collectLeadingXTerms(g_primeprime, ord) # g''_gamma * x^gamma, such that LT(g''_gamma) * x^gamma = LT(g'') 

        g_prime = g_prime + g_primeprime_gamma # Add the leading x-terms of g'' to g'
        g_primeprime = g_primeprime - g_primeprime_gamma # Subtract the leading x-terms of g'' from g''
    end

    #@req isInitiallyReduced(g_prime, ord) "g_prime=$g_prime is not initially reduced w.r.t. p-t under ord=$ord" # Check if g' is initially reduced w.r.t. p-t under ord
    return g_prime # Return g' which is initially reduced wrt p-t under > (no term of tail(g') is divisible by p) with <p-t, g'> = <p-t, g> and LT(g')= LT(g)
end

function sameDegreeReduce(G, p, ord)
    G = [ptReduce(g, p, ord) for g in G] # Apply ptReduce to each element of G
    G = sortByLeadingMonomial(G, ord) # Sort G by the leading term of each g w.r.t. ord

    @req length(G) > 0 "Number of terms is not greater than 0"
    @req !iszero(G[1]) "G[1] is zero"

    x_vars = gens(parent(G[1]))[2:end]
    # ----- FIRST PASS -----
    for i in 1:(length(G)-1)
        g_i = G[i] # Get the i-th element of G
        g_i_lt = lt(g_i, ord) # Get the leading term of g_i w.r.t. ord 

        exps = getExponents(g_i_lt) # Get the exponent vector of the leading term of g_i
        alpha_i = exps[2:end] # Get the X exponent vector of the leading term of g_i
        beta_i = exps[1] # Get the t exponent vector of the leading term of g_i

        x_alpha_i = prod(x_vars .^ alpha_i) # Calculate x^alpha
        Q, r = reduce_with_quotients(g_i, [x_alpha_i], ordering=ord) # Apply division algorithm to g_i w.r.t. x_vars

        @req length(Q) == 1 "Length of Q is not equal to 1"
        g_ialpha_i = Q[1] # Get the first element of Q

        for j in (i+1):length(G)
            g_j = G[j] # Get the j-th element of G
            Q, r = reduce_with_quotients(g_j, [x_alpha_i], ordering=ord) # Apply division algorithm to g_j w.r.t. g_i_lt

            @req length(Q) == 1 "Length of Q is not equal to 1"
            g_jalpha_i = Q[1]

            if !iszero(g_jalpha_i)
                g_j = (g_ialpha_i / t^(beta_i)) * g_j - (g_jalpha_i / t^(beta_i)) * g_i # Reduce g_j by g_i
                g_j = ptReduce(g_j, p, ord) # Reduce g_j w.r.t. p-t under >
                G[j] = g_j # Update G[j] with the initially reduced g_j
            end
        end
    end

    # ----- SECOND PASS -----
    for i in 1:(length(G)-1)
        g_i = G[i] # Get the i-th element of G

        for j in (i+1):length(G)
            g_j = G[j] # Get the j-th element of G
            g_j_lt = lt(g_j, ord) # Get the leading term of g_j w.r.t. ord 

            exps = getExponents(g_j_lt) # Get the exponent vector of the leading term of g_j
            alpha_j = exps[2:end] # Get the X exponent vector of the leading term of g_j
            beta_j = exps[1] # Get the t exponent vector of the leading term of g_j

            x_alpha_j = prod(x_vars .^ alpha_j) # Calculate x^alpha

            Q, r = reduce_with_quotients(g_j, [x_alpha_j], ordering=ord) # Apply division algorithm to g_j w.r.t. g_i_lt
            @req length(Q) == 1 "Length of Q is not equal to 1"
            g_jalpha_j = Q[1]

            Q, r = reduce_with_quotients(g_i, [x_alpha_j], ordering=ord) # Apply division algorithm to g_j w.r.t. g_i_lt
            @req length(Q) == 1 "Length of Q is not equal to 1"
            g_ialpha_j = Q[1]

            Q, r = reduce_with_quotients(g_ialpha_j, [t^(beta_j)], ordering=ord) # Apply division algorithm to g_jalpha_j w.r.t. g_i_lt
            if iszero(r) # Check if g_jalpha_j is divisible by t^beta_i
                #println("g_{$j}alpha_{$j} = ", g_jalpha_j, " is divisible by t^beta_{$j} = ", t^(beta_j))
                #println("g_{$i}alpha_{$j} = ", g_ialpha_j, " is divisible by t^beta_{$j} = ", t^(beta_j))
                g_i = (g_jalpha_j / t^(beta_j)) * g_i - (g_ialpha_j / t^(beta_j)) * g_j # Reduce g_j by g_i
                g_i = ptReduce(g_i, p, ord) # Reduce g_i w.r.t. p-t under >
                G[i] = g_i # Update G[i] with the initially reduced g_i
            end
        end
    end

    @req isInitiallyReduced(G, ord) "G=$G is not initially reduced w.r.t. G=$G, $p-t under ord=$ord" #(3) Check if G is initially reduced w.r.t. G, p-t under ord
    return G
end

function allAtOnceReduce(G, H, p, ord)
    # Apply Algorithm 4.5 from https://arxiv.org/abs/1512.02662 to reduce G w.r.t. H under ord
    # Input: G ∈ R[t, x1, ..., xn], H ∈ R[t, x1, ..., xn], p a prime number, ord a t-local monomial ordering on R
    # Output: H' ∈ R[t, x1, ..., xn] is H initially reduced w.r.t. G, H, p-t under ord

    d = sum(getXExponents(H[1])) # Get the degree of H[1] w.r.t. t
    for h in H
        @req lc(h, ord) == 1 "Leading coefficient of $h is not equal to 1 w.r.t. ord=$ord" #(3)
        for mono in monomials(h)
            @req sum(getXExponents(mono)) == d "Degree of $mono in $h is not equal to d=$d" #(1)
        end
    end
    for g in G
        for mono in monomials(g)
            @req sum(getXExponents(mono)) < d "Degree of $mono in $g is NOT < d=$d" #(2)
        end
    end
    H_lm = [lm(h, ord) for h in H] # Get the leading terms of H w.r.t. ord
    G_lm = [lm(g, ord) for g in G] # Get the leading terms of G w.r.t. ord
    @req length(H_lm) == length(unique(H_lm)) "H contains duplicate elements" #(4) Check if H contains duplicate elements
    @req intersect(H_lm, G_lm) == [] "H and G contain common elements" #(5) Check if H and G contain common elements

    H = sameDegreeReduce(H, p, ord) # Apply sameDegreeReduce to H w.r.t. p-t under ord
    E = []

    T = [] # Initialize T to an empty list
    for i in 1:length(H)
        h = H[i] # Get the i-th element of H
        h_lm = lm(h, ord) # Get the leading term of h w.r.t. ord
        h_collected = collectCoefficients(h, ord) # Collect the coefficients of h
        for (t_part, x_alpha) in h_collected
            h_alpha_lm = lm(t_part, ord) # Get the leading term of t_part of h_ialpha * x^alpha 
            if cmp(ord, h_alpha_lm * x_alpha, h_lm) == -1 # Check if LT(h_ialpha)x^alpha < LT(h_i)
                push!(T, (h_alpha_lm * x_alpha, i))
            end
        end
    end


    while !isempty(T)
        (s, i) = T[1]
        # Find the maximum (s,i) in T, with s maximal
        for j in 2:length(T)
            if cmp(ord, T[j][1], s) == 1
                (s, i) = T[j]
            end
        end

        for g in G
            g_lt = lt(g, ord) # Get the leading term of g w.r.t. ord
            Q, r = reduce_with_quotients(s, [g_lt], ordering=ord) # Apply division algorithm to g w.r.t. s
            if iszero(r) # Check if s is divisible by g
                s_lm = lm(s, ord) # Get the leading monomial of s w.r.t. ord
                g_lm = lm(g, ord) # Get the leading monomial of g w.r.t. ord
                push!(E, (s_lm / g_lm) * g) # Add the leading monomial of LM(s)/LM(g) * g to E
                H = sameDegreeReduce(union(E, H), p, ord) # Apply Algorithm 4.2 to E ∪ H
                T = filter(t -> cmp(t[1], s_lm) == -1, T) # Remove all elements of T with index i
            else
                T = filter(t -> t != (s, i), T) # Remove (s,i) from T
            end
        end
    end

    @req isInitiallyReduced(H, ord) "H=$H is not initially reduced w.r.t. G=$G, H=$H, $p-t under ord=$ord" #(3) Check if H is initially reduced w.r.t. G, H, p-t under ord
    return H # Return H which is initially reduced w.r.t. G, H, p-t under ord
end

function initialReduce(F, p, ord)
    # Apply Algorithm 4.7 from https://arxiv.org/abs/1512.02662 to reduce G w.r.t. H under ord
    # Input: F ∈ R[t, x1, ..., xn] x-homogeneous generating set, p a prime number, ord a t-local monomial ordering on R
    # Output: G ∈ R[t, x1, ..., xn] x-homogeneous & initially reduced groebner basis of <F> with respect to ord containing p-t
    for f in F
        d = sum(getXExponents(f)) # Get the degree of H[1] w.r.t. t
        for mono in monomials(f)
            @req sum(getXExponents(mono)) == d "Degree of $mono in $f is not equal to d=$d" #(1)
        end
    end

    I = ideal(F) # Create an ideal I from F
    G_primeprime = standard_basis(I, ordering=ord) # Get the standard basis of I w.r.t. ord
    for g_primeprime in G_primeprime
        d = sum(getXExponents(g_primeprime)) # Get the degree of H[1] w.r.t. t
        for mono in monomials(g_primeprime)
            @req sum(getXExponents(g_primeprime)) == d "Degree of $mono in $g_primeprime is not equal to d=$d" #(1)
        end
    end

    G_prime = [] # Initialize G' to an empty list
    for g in G_primeprime
        g_lc = lc(g, ord) # Get the leading term of g w.r.t. ord
        println("g_lc = ", g_lc) # Print the leading term of g w.r.t. ord
        if (g_lc != 1) # Check if g is not a unit and not divisible by p
            d, u, v = gcdx(g_lc, p) # Extended Euclid's algorithm to fing u,v such that u*g_lc + v*p = d
            g_prime = (u * g) + (v * lm(g, ord) * (p - t)) # Normalise g to have LC(g) = 1 
            push!(G_prime, g_prime) # Add g to G'
        end
    end

    if length(G_prime) > 1
        G_prime = minimise(G_prime, ord) # Minimise G' w.r.t. ord
    end

    G = []
    while !isempty(G_prime)
        #println("G_prime = ", G_prime) # Print G'
        d = minimum([sum(getXExponents(g)) for g in G_prime]) # Get the minimum degree of G'
        H_prime = filter(g -> sum(getXExponents(g)) == d, G_prime) # Get the elements of G' with degree d
        G_prime = filter(g -> sum(getXExponents(g)) > d, G_prime) # Remove the elements of G' with degree d
        H = allAtOnceReduce(G, H_prime, p, ord) # Apply Algorithm 4.5 to initially reduce H' w.r.t. G, H', p-t under ord
        G = union(G, H) # Add H to G
    end
    G = push!(G, p - t) # Add p - t to G
    return G # Return G which is initially reduced w.r.t. H under ord
end
function initially_reduce(F, p, ord)
    return initialReduce(F, p, ord) # Apply Algorithm 4.7 to reduce G w.r.t. H under ord
end