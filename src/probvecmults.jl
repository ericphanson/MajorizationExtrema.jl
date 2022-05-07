abstract type MajFlowNorm end

struct OneNorm <: MajFlowNorm end

struct InfNorm <: MajFlowNorm end

(nrm::InfNorm)(s, q) = norm(s - q, Inf)

export SortedProbVecMult

Base.:(-)(spvm::SortedProbVecMult, spvm2::SortedProbVecMult) = SortedProbVecMult(collect(spvm) - collect(spvm2))
struct SortedProbVecMult{T, VT <: AbstractVector{T}} <: AbstractVector{T}
    distinct_entries::VT
    multiplicities::Vector{Int}
end


Base.size(r::SortedProbVecMult) = (sum(r.multiplicities),)
Base.eltype(r::SortedProbVecMult{T}) where {T} = T

Base.@propagate_inbounds function Base.getindex(r::SortedProbVecMult, i::Int)
    @boundscheck (1 <= i <= length(r)) || throw(BoundsError(r,i))
    ks = r.multiplicities
    seen_so_far = 0
    for j = eachindex(ks)
        seen_so_far += ks[j]
        if seen_so_far >= i
            return r.distinct_entries[j]
        end
    end
end

function SortedProbVecMult(v::AbstractVector)
    v_sorted = sort(v; rev = true)
    distinct_entries = eltype(v)[]
    multiplicities = Int[]
    for elt in v_sorted
        if isempty(distinct_entries) || distinct_entries[end] != elt
            push!(distinct_entries, elt)
            push!(multiplicities, 1)
        else
            multiplicities[end] += 1
        end

    end
    SortedProbVecMult{eltype(distinct_entries), typeof(distinct_entries)}(distinct_entries, multiplicities)
end


function Base.collect(spvm::SortedProbVecMult)
    @unpack distinct_entries, multiplicities = spvm
    v = eltype(distinct_entries)[]
    sizehint!(v, length(spvm))
    for (i, elt) in enumerate(distinct_entries)
        for j = 1:multiplicities[i]
            push!(v, elt)
        end
    end
    return v
end

# Default
majmin(v::SortedProbVecMult, ϵ) = majmin(OneNorm(), v, ϵ)

function majmin(nrm::MajFlowNorm, spvm::SortedProbVecMult, ϵ)
    majmin!(nrm, deepcopy(spvm), ϵ)
end


function majmin!(nrm::OneNorm, spvm::SortedProbVecMult, ϵ)
    @unpack distinct_entries, multiplicities = spvm
    @assert ϵ >= 0
    if length(distinct_entries) == 1
        return spvm
    elseif length(distinct_entries) == 2
        _majmin_2_entries(spvm::SortedProbVecMult, ϵ)
    else
        kp = first(multiplicities)
        km = last(multiplicities)

        δp = kp * (distinct_entries[1] - distinct_entries[2])
        δm = km * (distinct_entries[end-1] - distinct_entries[end])
        gap = min(δp, δm)
        if ϵ < gap
            distinct_entries[1] -= ϵ / kp
            distinct_entries[end] += ϵ / km
            return spvm
        else
            if δp < δm
                # evolve for time δp and recurse
                deleteat!(distinct_entries, 1)
                deleteat!(multiplicities, 1)
                multiplicities[1] += kp
                distinct_entries[end] += δp / km
                return majmin!(nrm, spvm, ϵ - δp)
            elseif δm < δp
                # evolve for time δm and recurse
                pop!(distinct_entries)
                pop!(multiplicities)
                multiplicities[end] += km
                distinct_entries[1] -= δm / kp
                return majmin!(nrm, spvm, ϵ - δm)
            else
                # evolve for time gap=δm=δp and recurse
                deleteat!(distinct_entries, 1)
                deleteat!(multiplicities, 1)
                multiplicities[1] += kp
                pop!(distinct_entries)
                pop!(multiplicities)
                multiplicities[end] += km
                return majmin!(nrm, spvm, ϵ - gap)
            end
        end
    end
end


# This actually works for either 1-norm or infinity norm (I think)
# With only 2 entries, the next crossing isn't just determined by
# gap = min(δp, δm).
function _majmin_2_entries(spvm::SortedProbVecMult, ϵ)
    @unpack distinct_entries, multiplicities = spvm
    @assert length(distinct_entries) == length(multiplicities) == 2
    kp, km = multiplicities
    μp, μm = distinct_entries
    # the larger entry is decreasing at a rate t/kp
    # the smaller entry is increasing at a rate t/km
    # Need to solve   `μp - t/kp == μm + t/km` for `t`
    # t = (μp - μm) / ( 1/km + 1/kp )
    t = (μp - μm) * (km * kp) / (km + kp)
    if ϵ < t
        distinct_entries[1] -= ϵ / kp
        distinct_entries[2] += ϵ / km
        return spvm
    else
        pop!(distinct_entries)
        pop!(multiplicities)
        distinct_entries[] -= t / kp
        multiplicities[] += km
        return spvm
    end

end

function majmin!(nrm::InfNorm, spvm::SortedProbVecMult, ϵ)
    @unpack distinct_entries, multiplicities = spvm
    @assert issorted(distinct_entries; rev=true)
    @assert sum(spvm) == 1
    T = eltype(spvm)
    @assert T == Rational{BigInt} == typeof(ϵ)
    # ϵ = min(ϵ, nrm(spvm, ones(T, length(spvm)) .// length(spvm)))
    m = length(distinct_entries)
    # @assert m > 2
    if m == 2
        @show Float64.(distinct_entries), Float64(ϵ), multiplicities
        return _majmin_2_entries(spvm, ϵ)
    end

    # Algorithm for computing the majorization minimizer over an infinity-norm ball on the simplex:
    # 1. The largest half of distinct entries decrease infinitesimally with a rate proportional to their multiplicity
    # 2. The smallest half of distinct entries increase infinitesimally with a rate proportional to their multiplicity
    # 3. If there's an odd number of distinct entries, the middle one doesn't move
    # The number of distinct entries is "re-evaluated infinitesimally"
    # in reality we just compute the crossings and re-evaluate then.

    # Compute the gaps
    # the larger entry is decreasing at a rate t/kp
    # the smaller entry is increasing at a rate t/km
    # Need to solve   `μp - t/kp == μm + t/km` for `t`
    # t = (μp - μm) / ( 1/km + 1/kp )
    gap(i, j) = (distinct_entries[i] - distinct_entries[j]) * (multiplicities[i] * multiplicities[j]) / (multiplicities[i] + multiplicities[j])

    # the larger entry is decreasing at a rate t/kp
    # the smaller entry is not moving
    # Need to solve  `μp - t/kp == μm` for `t`
    # t = (μp - μm)*kp
    gap2(i, j) = (distinct_entries[i] - distinct_entries[j]) * multiplicities[i]


    # the larger entry is not moving
    # the smaller entry is increasing at rate t/km
    # Need to solve  `μp  == μm + t/km` for `t`
    # t = (μp - μm)*km
    gap3(i, j) = (distinct_entries[i] - distinct_entries[j]) * multiplicities[j]

    if isodd(m)
        k = (m-1)÷2
        X = [-ones(T, k); zero(T); ones(T, k)]
        @assert length(X) == m
        @assert sum(X) == 0
        X = X .// multiplicities
        ϵ₁ = min(gap2(k, k+1), gap3(k+1, k+2), gap(k, k+2))
        @assert ϵ₁ >= 0
    else
        k = m÷2
        X = [-ones(T, k); ones(T, k)]
        @assert length(X) == m
        @assert sum(X) == 0
        X = X .// multiplicities
        ϵ₁ = gap(k, k+1)
        @assert ϵ₁ >= 0
    end
    @show Float64.(distinct_entries), Float64(ϵ), Float64(ϵ₁), multiplicities

    @assert !iszero(ϵ₁)
    # iszero(ϵ₁) && return spvm
    if ϵ₁ < ϵ
        remaining_ϵ = ϵ - ϵ₁
        @. spvm.distinct_entries += X*ϵ₁
        @assert issorted(distinct_entries; rev=true)

        spvm = SortedProbVecMult(collect(spvm))
        return majmin!(nrm, spvm, remaining_ϵ)
    else
        @. spvm.distinct_entries += X*ϵ
        @assert issorted(distinct_entries; rev=true)

        spvm = SortedProbVecMult(collect(spvm))
        
        @show Float64.(distinct_entries), Float64(ϵ), multiplicities
        return spvm
    end

end


using Convex, Tulip
function majmin_inf(p, ϵ)
    ϵ < 1e-12 && return p
    perm = sortperm(p)
    p = p[perm]
    d = length(p)
    q = Variable(d)
    constraints = [ q >= 0, sum(q) == 1, - ϵ <= q - p, q - p <= ϵ]
    y = [big(0.0)]
    for k = 1:d
        prob = minimize(sumlargest(q, k), constraints; numeric_type=BigFloat)
        solve!(prob, Tulip.Optimizer{BigFloat}(); silent_solver=true)
        push!(y, prob.optval)
    end
    return sort(diff(y))[invperm(perm)]
end

using Test
function test_inf(v, ϵ)
    sort!(v; rev=true)
    v = Rational{BigInt}.(v)
    ϵ = Rational{BigInt}(ϵ)
    @assert sum(v) == 1
    @assert 0 < ϵ < 1

    x = majmin(InfNorm(), SortedProbVecMult(v), ϵ)
    y = MajorizationExtrema.majmin_inf(v, ϵ)
    @test Float64.(x) ≈ Float64.(y)
end

# TODO:
# 1. Fix m=2 case
# 2. Prove formula/algorithm
# 3. Prove semigroup property
