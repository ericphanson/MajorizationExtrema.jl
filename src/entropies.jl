function entropy(q::AbstractVector)
    eta(x) = x <= 0 ? zero(x) : -x*log2(x)
    return sum(eta, q)
end

entropy(ρ::AbstractMatrix) = entropy(eigvals(ρ))