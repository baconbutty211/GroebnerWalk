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
                println("g_{$j}alpha_{$j} = ", g_jalpha_j, " is divisible by t^beta_{$j} = ", t^(beta_j))
                println("g_{$i}alpha_{$j} = ", g_ialpha_j, " is divisible by t^beta_{$j} = ", t^(beta_j))
                g_i = (g_jalpha_j / t^(beta_j)) * g_i - (g_ialpha_j / t^(beta_j)) * g_j # Reduce g_j by g_i
                g_i = ptReduce(g_i, p, ord) # Reduce g_i w.r.t. p-t under >
                G[i] = g_i # Update G[i] with the initially reduced g_i
            end
        end
    end

    return G
end