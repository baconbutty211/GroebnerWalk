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

function getXExponents(term::MPolyRingElem)
    # Get the exponents of x1, ..., xn from a term in R[t, x1, ..., xn]
    # Input: term ∈ R[t, x1, ..., xn]
    # Output: exponents ∈ ZZ^n is the exponent vector of x1, ..., xn
    @req length(exponents(term)) > 0 "Length of exponents is not greater than 0"
    return collect(exponents(t))[1][2:end] # Return the exponent vector excluding t
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