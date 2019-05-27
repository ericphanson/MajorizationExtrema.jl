"""
    TV(p, q)

Returns the total variation distance between `p` and `q` (half the sum of the absolute value of the differences).
"""
TV(p::AbstractVector, q::AbstractVector) = (1/2)*norm(p - q, 1)

function TV(p::AbstractVector{T1}, q::AbstractVector{T2}) where {T1 <: Rational, T2 <: Rational}
    length(p) == length(q) || throw(DimensionMismatch("Input vectors must be of the same length"))
    return (1//2)*sum(abs, p - q)
end

"""
    ≺(p::AbstractVector{T1}, q::AbstractVector{T2}; tol = (T1 <: AbstractFloat || T2 <: AbstractFloat) ? T2(1e-8) : zero(T2) ) where {T1, T2}

Returns true if `q` majorizes `p` and false otherwise. The keyword argument `tol` specifies a tolerance for the comparisons. Can be used in infix form, i.e. `p ≺ q`.
"""
function ≺(p::AbstractVector{T1}, q::AbstractVector{T2}; tol = (T1 <: AbstractFloat || T2 <: AbstractFloat) ? T2(1e-8) : zero(T2) ) where {T1, T2}
    length(p) == length(q) || throw(DimensionMismatch("Input vectors must be of the same length"))
    pv = sort(p, rev = true)
    qv = sort(q, rev = true)
    isapprox(sum(p), sum(q), atol=tol, rtol=0) || return false
    return all(cumsum(pv) .<= ( cumsum(qv) .+ tol ) )
end


"""
    majmax(q::AbstractVector, ϵ)

Returns the maximum in majorization order over the total variation ball (of probability vectors) of radius `ϵ` around a probability vector `q`.
"""
function majmax(q::AbstractVector, ϵ)
    d = length(q)
    T = eltype(q[1] + ϵ)
    inds = sortperm(q)
    q = q[inds]
    p = zeros(T,d)
    budget = ϵ

    p[1] = max(zero(T), q[1] - budget )
    budget = budget - abs(q[1]-p[1])

    for j = 2:d-1
        p[j] = max(zero(T), q[j] - budget)
        budget = budget - abs(q[j]-p[j])
    end
    p[d] = q[d] + ϵ - budget
    invpermute!(p, inds)
    return p
end

@inline function get_alpha1(qup::AbstractVector, ϵ)
    d = length(qup)
    cs = qup[1] + ϵ
    for j = 1:d-1
        if cs <= j*qup[j+1]
            return j, cs/j
        end
        cs += qup[j+1]
    end
    return nothing
end

@inline function get_alpha2(qup::AbstractVector, ϵ)
    d = length(qup)
    cs = qup[d] - ϵ
    for (idx, j) = enumerate(d:-1:2)
        if cs >= idx*qup[j-1]
            return idx, cs/idx
        end
        cs += qup[j-1]
    end
    return nothing
end

uniform(T, d) = inv(T(d)) * ones(T, d)

"""
    majmin(q::AbstractVector, ϵ)

Returns the minimum in majorization order over the total variation ball (of probability vectors) of radius `ϵ` around a probability vector `q`.
"""
function majmin(q::AbstractVector, ϵ)
    T = typeof((q[1] + ϵ)/2)
    d = length(q)
    inds = sortperm(q, rev = false)

    qup = convert.(T, q[inds])
    VT = typeof(qup)

    n1a1 = get_alpha1(qup, ϵ)
    n2a2 = get_alpha2(qup, ϵ)

    if n1a1 === nothing || n2a2 === nothing        
         return VT(uniform(T, d))
    end

    n1, alpha1 = n1a1
    n2, alpha2 = n2a2

    if alpha1 > 1/d || alpha2 < 1/d
        return VT(uniform(T,d))
    end

    out = qup
    out[1:n1] .= alpha1
    out[d-n2+1:d] .= alpha2
    invpermute!(out, inds)
    return VT(out)
end

