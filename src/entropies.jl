function entropy(q::AbstractVector)
    η = x -> x <= 0 ? zero(x) : -x * log2(x)
    return sum(η, q)
end

entropy(ρ::AbstractMatrix) = entropy(eigvals(ρ))

renyi_entropy(α::Real) = r -> renyi_entropy(r, α)

function renyi_entropy(q::AbstractVector, α)
    α == 1 && return entropy(q)
    α == Inf && return -log2(maximum(q))
    α == 0 && return sum(q .> 0)
    φ = x -> x <= 0 ? zero(x) : x^α
    return inv(1 - α) * log2(sum(φ, q))
end

tsallis_entropy(q::Real) = r -> tsallis_entropy(r, q)

function tsallis_entropy(r::AbstractVector, q)
    q == 1 && return entropy(r)
    φ = x -> x <= 0 ? zero(x) : x^q
    return (sum(φ, r) - 1) * inv(1 - q)
end

subentropy(r::AbstractVector) = subentropy(SortedProbVecMult(r))

function subentropy(spvm::SortedProbVecMult)
    μs = spvm.distinct_entries
    ks = spvm.multiplicities
    n = sum(ks)
    s = zero(eltype(μs))
    for i in eachindex(μs)
        iszero(μs[i]) && continue
        s -= _subentropy_fi(spvm, i, ks[i] - 1, μs[i]) / factorial(ks[i] - 1)
    end
    return s
end

# this is slow for large `n`
# computes the `n`th derivative of `f` at `x`
function _deriv(f, n, x)
    if n == 0
        return f(x)
    end
    D = f -> (s -> ForwardDiff.derivative(f, s))
    Dn = foldl(∘, (D for _ = 1:n))
    Dn(f)(x)
end

function _subentropy_fi(spvm, i, k, x)
    μs = spvm.distinct_entries
    ks = spvm.multiplicities
    d = sum(ks)
    fi = x -> begin
        itr = ((x - μs[j])^ks[j] for j = eachindex(μs) if j ≠ i)
        prod_others = isempty(itr) ? one(x) : prod(itr)
        x^d * log(x) / prod_others
    end
    _deriv(fi, k, x)
end
