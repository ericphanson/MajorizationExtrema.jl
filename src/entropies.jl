function entropy(q::AbstractVector)
    η = x -> x <= 0 ? zero(x) : -x*log2(x)
    return sum(η, q)
end

entropy(ρ::AbstractMatrix) = entropy(eigvals(ρ))

renyi_entropy(α::Real) = r -> renyi_entropy(r, α)

function renyi_entropy(q::AbstractVector, α)
    α == 1 && return entropy(q)
    α == Inf && return - log2(maximum(q))
    α == 0 && return sum(q .> 0)
    φ = x -> x <= 0 ? zero(x) : x^α
    return inv(1-α) * log2( sum(φ, q))
end

tsallis_entropy(q::Real) = r -> tsallis_entropy(r, q)

function tsallis_entropy(r::AbstractVector, q)
    q == 1 && return entropy(r)
    φ = x -> x <= 0 ? zero(x) : x^q
    return (sum(φ, r) - 1)*inv(1-q)
end
