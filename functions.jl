using Oscar

function lt(f::MPolyRingElem, ord::MonomialOrdering)
    # Get the leading term of a polynomial f w.r.t. an ordering ord
    # Input: f ∈ R[t, x1, ..., xn], ord a t-local monomial ordering on R
    # Output: lt(f) ∈ R[t, x1, ..., xn] is the leading term of f w.r.t. ord
    @req length(f) > 0 "Number of terms is not greater than 0"

    lead_term = term(f, 1) # Get the first monomial of f
    lead_mono = monomial(f, 1) # Get the first monomial of f
    for i in 2:length(f)
        t = term(f, i) # Extract the monomial part of the term
        m = monomial(f, i) # Extract the monomial part of the term
        if cmp(ord, m, lead_mono) == 1
            lead_mono = m # Update the leading monomial
            lead_term = t # Update the leading term
        end
    end
    return lead_term
end
function lm(f::MPolyRingElem, ord::MonomialOrdering)
    monomial(lt(f, ord), 1)
end
function lc(f::MPolyRingElem, ord::MonomialOrdering)
    coeff(lt(f, ord), 1) # Get the coefficient of the leading term of f w.r.t. ord
end

function getExponents(term::MPolyRingElem)
    # Get the exponents of x1, ..., xn from a term in R[t, x1, ..., xn]
    # Input: term ∈ R[t, x1, ..., xn]
    # Output: exponents ∈ ZZ^n is the exponent vector of x1, ..., xn
    @req length(exponents(term)) > 0 "Length of exponents is not greater than 0"
    return collect(exponents(term))[1] # Return the exponent vector excluding t
end
function getXExponents(term::MPolyRingElem)
    # Get the exponents of x1, ..., xn from a term in R[t, x1, ..., xn]
    # Input: term ∈ R[t, x1, ..., xn]
    # Output: exponents ∈ ZZ^n is the exponent vector of x1, ..., xn
    return getExponents(term)[2:end] # Return the exponent vector excluding t
end
function getTExponent(term::MPolyRingElem)
    # Get the exponents of x1, ..., xn from a term in R[t, x1, ..., xn]
    # Input: term ∈ R[t, x1, ..., xn]
    # Output: exponents ∈ ZZ^n is the exponent vector of x1, ..., xn
    return getExponents(term)[1] # Return the exponent vector excluding t
end

function getTCoefficient(term::MPolyRingElem)
    # Get the coefficient from a term in R[t, x1, ..., xn]
    # Input: term ∈ R[t, x1, ..., xn]
    # Output: coeff ∈ R[t] is the coefficient of term
    @req length(exponents(term)) > 0 "Length of exponents is not greater than 0"
    c = coeff(term, 1) # Get the coefficient of the term
    @req c !== 0 "Coefficient is zero"

    alpha_prime = collect(exponents(t))[1][1] # Return the coefficient of t
    return c * t^alpha_prime # Return the coefficient in R[t]
end

function collectXTerms(f::MPolyRingElem, alpha::Vector{Int})
    # Assume t variable is the first variable in the polynomial ring
    # Collect the terms of f with same  x1, ..., xn exponent vector alpha
    # Input: f ∈ R[t, x1, ..., xn]
    # Output: f_alpha ∈ R[t][x1,...,xn] is the monomial in R[t][x1,...,xn] with exponent vector alpha
    @req length(f) > 0 "Number of terms is not greater than 0"
    @req length(alpha) == nvars(parent(f)) - 1 "Length of alpha is not equal to 1 less than the number of variables in the polynomial ring"
    g = 0 # Initialize g to 0
    for t in terms(f)
        if collect(exponents(t))[1][2:end] == alpha
            g += t
        end
    end
    return g # Return the sum of the terms with same exponent vector alpha
end
function collectLeadingXTerms(f::MPolyRingElem, ord::MonomialOrdering)
    # Assume t variable is the first variable in the polynomial ring
    # Collect the terms of f with same  x1, ..., xn exponent vector alpha
    # Input: f ∈ R[t, x1, ..., xn]
    # Output: f_alpha ∈ R[t][x1,...,xn] is the monomial in R[t][x1,...,xn] with exponent vector alpha
    @req length(f) > 0 "Number of terms is not greater than 0"
    alpha = getXExponents(lt(f, ord)) # Get the exponent vector of the leading term of f w.r.t. lex ordering
    return collectXTerms(f, alpha) # Return the sum of the terms with same exponent vector alpha
end

function collectLeadingXCoefficients(f::MPolyRingElem, ord::MonomialOrdering)
    # Assume t variable is the first variable in the polynomial ring
    # Collect the leading coefficient of f w.r.t. ord
    # Input: f ∈ R[t, x1, ..., xn]
    # Output: lc_alpha ∈ R[t] is the leading coefficient of f w.r.t. ord
    @req length(f) > 0 "Number of terms is not greater than 0"
    g = collectLeadingXTerms(f, ord) # Collect the leading x-terms of f w.r.t. ord
    alpha = getXExponents(lt(g, ord)) # Get the exponent vector of the leading term of g w.r.t. lex ordering
    x_vars = gens(parent(g))[2:end] # Get the variables of the polynomial ring
    x_alpha = prod(x_vars .^ alpha)
    Q, r = reduce_with_quotients(g, x_alpha, ordering=ord) # Apply division algorithm to g w.r.t. x_vars to get {q1, ..., qn} (in finite field for performance reasons)

    @req iszero(r) "Remainder is not zero"
    @req length(Q) == 1 "Length of Q is not equal to 1"

    return Q[1]
end

function collectCoefficients(f::MPolyRingElem, ord::MonomialOrdering)
    # Assume t variable is the first variable in the polynomial ring
    # Input: f ∈ R[t, x1, ..., xn]
    # Output: g ∈ R[t][x1, ..., xn]
    @req length(f) > 0 "Number of terms is not greater than 0"
    f_prime = [] # Initialize f' to lt(f, ord)
    f_primeprime = f # Initialize f'' to f - lt(f, ord)

    x_vars = gens(parent(f))[2:end] # Get the variables of the polynomial ring
    while !iszero(f_primeprime)
        f_lt = lt(f_primeprime, ord) # Get the leading term of f w.r.t. ord
        alpha = getXExponents(f_lt) # Get the exponent vector of the leading term of term w.r.t. lex ordering
        x_alpha = prod(x_vars .^ alpha) # Calculate x^alpha
        Q, r = reduce_with_quotients(f_primeprime, [x_alpha], ordering=ord) # Apply division algorithm to g w.r.t. x_vars to get {q1, ..., qn} (in finite field for performance reasons)

        @req length(Q) == 1 "Length of Q is not equal to 1"

        pushfirst!(f_prime, (Q[1], x_alpha)) # Add (f_alpha, x^alpha) to f'
        f_primeprime = r # Subtract the remainder from f''
    end
    return f_prime # Return f' which is the list of the terms f_alpha(t) * x^alpha in f in the form [(f_alpha, x^alpha), ...] 
end

function sortByLeadingMonomial(F, ord::MonomialOrdering, ascending::Bool=false)
    # Sort the elements of F by their leading monomial w.r.t. ord
    # Input: F ∈ R[t, x1, ..., xn]
    # Output: F_sorted ∈ R[t, x1, ..., xn] is the sorted list of elements of F by their leading monomial w.r.t. ord
    @req length(F) > 0 "Number of terms is not greater than 0"
    @req length(F) == length(unique(F)) "F contains duplicate elements"

    F_sorted = []
    for f in F
        if length(F_sorted) == 0
            push!(F_sorted, f) # Add the first element to the sorted list
        else
            i = 1 # Initialize index to 1
            if (ascending)
                while i <= length(F_sorted) && cmp(ord, lm(f, ord), lm(F_sorted[i], ord)) == 1
                    i += 1 # Increment index until the correct position is found
                end
            else
                while i <= length(F_sorted) && cmp(ord, lm(f, ord), lm(F_sorted[i], ord)) == -1
                    i += 1 # Increment index until the correct position is found
                end
            end
            insert!(F_sorted, i, f) # Insert f at the correct position in the sorted list
        end
    end
    return F_sorted # Return the sorted list of elements of F by their leading monomial w.r.t. ord
end

function minimise(G, ord)
    # Minimise G w.r.t. ord
    # Input: G ∈ R[t, x1, ..., xn], ord a t-local monomial ordering on R
    # Output: G_min ∈ R[t, x1, ..., xn] is the minimised list of elements of G w.r.t. ord
    @req length(G) > 0 "Number of terms is not greater than 0"
    G = sortByLeadingMonomial(G, ord, true) # Sort G by their leading monomial w.r.t. ord
    G_min = [] # Initialize G_min to empty list

    G_lm = [(lm(G[i], ord), g) for g in G] # Get the leading terms of G w.r.t. ord
    for (g_lm, i) in G_lm
        G_minus_g_lm = filter(_lm -> _lm != g_lm, G_lm) # Get the leading terms of G without g
        Q, r = reduce_with_quotients(g_lm, G_minus_g_lm, ordering=ord)
        if r == g_lm # Check if g divides a term in g'
            push!(G_min, g) # Add g to G_min if it does not divide a term in g'
        end
    end
    return G_min # Return the minimised list of elements of G w.r.t. ord
end

function isMinimised(G, ord)
    # Minimise G w.r.t. ord
    # Input: G ∈ R[t, x1, ..., xn], ord a t-local monomial ordering on R
    # Output: G_min ∈ R[t, x1, ..., xn] is the minimised list of elements of G w.r.t. ord
    @req length(G) > 0 "Number of terms is not greater than 0"
    G = sortByLeadingMonomial(G, ord, true) # Sort G by their leading monomial w.r.t. ord

    G_lm = [lm(G[i], ord) for g in G] # Get the leading terms of G w.r.t. ord
    for g_lm in G_lm
        G_minus_g_lm = filter(_lm -> _lm != g_lm, G_lm) # Get the leading terms of G without g
        Q, r = reduce_with_quotients(g_lm, G_minus_g_lm, ordering=ord)
        if r != g_lm
            return false
        end
    end
    return true # Return the minimised list of elements of G w.r.t. ord
end
function isReduced(G, ord)
    # Check if G is reduced w.r.t. ord
    # Input: G ∈ R[t, x1, ..., xn], ord a t-local monomial ordering on R
    # Output: true if G is reduced w.r.t. ord, false otherwise
    @req length(G) > 0 "Number of terms is not greater than 0"
    G_lm = [lm(_g, ord) for _g in G] # Get the leading terms of G w.r.t.
    for g in G
        g_lm = lm(g, ord) # Get the leading term of g w.r.t. ord
        G_minus_g_lm = filter(_lm -> _lm != g_lm, G_lm) # Get the leading terms of G without g
        Q, r = reduce_with_quotients(g, G_minus_g_lm, ordering=ord)
        #println("g = ", g, ", G_minus_g_lm = ", G_minus_g_lm, ", Q = ", Q, ", r = ", r)
        if r != g # Check if g divides a term in g'
            return false
        end
    end
    return true # G is reduced w.r.t. ord
end
function isInitiallyReduced(G, ord)
    # Check if G is initially reduced w.r.t. p-t, H, G under >
    # Input: G ∈ R[t, x1, ..., xn], H ∈ R[t, x1, ..., xn], p a prime number, > a t-local monomial ordering on R
    # Output: true if G is initially reduced w.r.t. p-t under >, false otherwise
    @req length(G) > 0 "Number of terms is not greater than 0"
    @req length(H) > 0 "Number of terms is not greater than 0"

    G_prime = []
    for g in G
        g_collected = collectCoefficients(g, ord) # Collect the coefficients of g w.r.t. ord
        g_prime = 0
        for (g_alpha, x_alpha) in g_collected
            g_alpha_lt = lt(g_alpha, ord) # Get the leading term of g_alpha w.r.t. ord
            g_prime += g_alpha_lt * x_alpha # Add the coefficients of g w.r.t. ord
        end
        push!(G_prime, g_prime) # Add the leading term of g to G'
    end

    return isReduced(G_prime, ord) # G is initially reduced w.r.t. p-t under >
end

function isXHomogeneous(f)
    # Check if f is x-homogeneous
    # Input: f ∈ R[t, x1, ..., xn]
    # Output: true if f is x-homogeneous, false otherwise

    d = sum(getXExponents(f)) # Get the degree of the first monomial of f
    for mono in monomials(f)
        if sum(getXExponents(mono)) != d
            return false
        end
    end
    return true # f is x-homogeneous
end
function isXHomogeneous(F)
    # Check if f is x-homogeneous
    # Input: f ∈ R[t, x1, ..., xn]
    # Output: true if f is x-homogeneous, false otherwise

    d = sum(getXExponents(F[1])) # Get the degree of the first monomial of f
    for f in F
        for mono in monomials(f)
            if sum(getXExponents(mono)) != d
                return false
            end
        end
    end
    return true # f is x-homogeneous
end