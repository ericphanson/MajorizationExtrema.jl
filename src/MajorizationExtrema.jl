module MajorizationExtrema

using LinearAlgebra, Random

export simplexpt, randsimplexpt
export TV, ≺
export majmax, majmin, localbound

include("fisher_yates!.jl")

"""
    simplexpt(unif)

Takes a vector of length `d-1` of numbers between `0` and `1` and converts it a point on the standard `d` dimensional simplex.
"""
function simplexpt(unif)
    d = length(unif) + 1; T = eltype(unif)
    w= zeros(T,d+1)
    w[2:d] .= sort(unif)
    w[d+1] = one(T)
    diff(w)
end

"""
    randsimplexpt(d)

Generates points uniformly at random on the standard `d-1` dimensional simplex using an algorithm by [Smith and Tromble](http://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf).
"""
randsimplexpt(d)  =  simplexpt(rand(Float64, d-1))

"""
    randsimplexpt(d, N::T) where {T <: Integer}

Generates rational numbers with denominator at most `N` uniformly at random on the `d-1` dimensional simplex using an algorithm by [Smith and Tromble](http://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf).
"""
function randsimplexpt(d, N::T) where {T <: Integer}
    unif = zeros(T,d-1)
    fisher_yates_sample!(1:N,unif)
    unif = unif .// N
    simplexpt(unif)
end

function TV(p::AbstractVector{T1}, q::AbstractVector{T2}) where {T1 <: Rational, T2 <: Rational}
    length(p) == length(q) || throw(ArgumentError("Input vectors must be of the same length"))
    return (1//2)*sum(abs(p[j] - q[j]) for j = eachindex(p))
end

"""
    TV(p, q)

Returns the total variation distance between `p` and `q` (half the sum of the absolute value of the differences).
"""
TV(p, q) = (1/2)*norm(p - q, 1)


"""
    ≺(p::AbstractVector{T1}, q::AbstractVector{T2}; tol = (T1 <: AbstractFloat || T2 <: AbstractFloat) ? T2(1e-8) : zero(T2) ) where {T1, T2}

Returns true if `q` majorizes `p` and false otherwise. The keyword argument `tol` specifies a tolerance for the comparisons. Can be used in infix form, i.e. `p ≺ q`.
"""
function ≺(p::AbstractVector{T1}, q::AbstractVector{T2}; tol = (T1 <: AbstractFloat || T2 <: AbstractFloat) ? T2(1e-8) : zero(T2) ) where {T1, T2}
    length(p) == length(q) || throw(ArgumentError("Input vectors must be of the same length"))
    pv = sort(p, rev = true);
    qv = sort(q, rev = true);
    isapprox(sum(p), sum(q), atol=tol, rtol=0) || return false
    return all(cumsum(pv) .<= ( cumsum(qv) .+ tol ) )
end


"""
    majmax(q, ϵ)

Returns the maximum in majorization order over the total variation ball (of probability vectors) of radius `ϵ` around a probability vector `q`.
"""
function majmax(q, ϵ)
    d = length(q);
    T = eltype(q);
    inds =sortperm(q);
    q = q[inds];
    p = zeros(T,d);
    budget = ϵ;

    p[1] = max(0, q[1] - budget );
    budget = budget - abs(q[1]-p[1]);

    for j = 2:d-1
        p[j] = max(0, q[j] - budget);
        budget = budget - abs(q[j]-p[j]);
    end
    p[d] = q[d] + ϵ - budget;

    return p[invperm(inds)]
end

@inline function make_alphas(q, ϵ)
    d = length(q)
    return (cumsum(q) .+ ϵ) ./ (1:d)
end

"""
    majmin(q, ϵ)

Returns the minimum in majorization order over the total variation ball (of probability vectors) of radius `ϵ` around a probability vector `q`.
"""
function majmin(q, ϵ::ϵT) where {ϵT}
    T = eltype(q)
    d = length(q)
    inds =sortperm(q, rev = false);

    qup = q[inds]

    alpha1s = make_alphas(qup, ϵ);

    n1 =  findfirst(alpha1s[1 : d-1]  .<= qup[2 : d]);
    n1 === nothing && return inv(T(d)) * ones(T,d);

    alpha1 = alpha1s[n1];
    alpha1 > 1/d && return inv(T(d)) * ones(T,d);

    qdown = sort(q, rev=true);
    alpha2s = make_alphas(qdown, -1*ϵ);

    n2 = findfirst(alpha2s[1 : d-1] .>= qdown[2 : d] );
    n2 === nothing && return inv(T(d)) * ones(T,d);

    alpha2 = alpha2s[n2];
    alpha2 < 1/d && return inv(T(d)) * ones(T,d);

    out = qup
    out[1:n1] .= alpha1
    out[d-n2+1:d] .= alpha2

    return out[invperm(inds)]
end

"""
    localbound(f, p, ϵ)

Returns `δ` such that `|f(p) - f(q)| <= δ` for any `q` with `TV(p, q) <= ϵ`, where `p` and `q` are probability vectors. Requires `f` to be Schur convex or Schur concave (but does not check this condition).
"""
function localbound(f, p, ϵ)
    fp = f(p)
    f_max = f(majmax(p, ϵ))
    f_min = f(majmin(p, ϵ))
    return max(abs(fp - f_max), abs(fp - f_min))
end

end # module
