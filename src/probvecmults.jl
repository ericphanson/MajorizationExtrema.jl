abstract type MajFlowNorm end

struct OneNorm <: MajFlowNorm end

struct InfNorm <: MajFlowNorm end

export SortedProbVecMult
struct SortedProbVecMult{VT}
    distinct_entries::VT
    multiplicities::Vector{Int}
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
    SortedProbVecMult{typeof(distinct_entries)}(distinct_entries, multiplicities)
end


function Base.collect(spvm::SortedProbVecMult)
    @unpack distinct_entries, multiplicities = spvm
    v = eltype(distinct_entries)[]

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
       _majmin_onenorm_2entries!(spvm::SortedProbVecMult, ϵ)
    else
        kp = first(multiplicities)
        km = last(multiplicities)

        δp = kp*(distinct_entries[1] - distinct_entries[2])
        δm = km*(distinct_entries[end-1] - distinct_entries[end])
        gap = min(δp, δm)
        if ϵ < gap
            distinct_entries[1] -= ϵ/kp
            distinct_entries[end] += ϵ/km
            return spvm
        else
            if δp < δm
                # evolve for time δp and recurse
                deleteat!(distinct_entries, 1)
                deleteat!(multiplicities, 1)
                multiplicities[1] += kp
                distinct_entries[end] += δp/km
                return majmin!(nrm, spvm, ϵ - δp)
            elseif δm < δp
                # evolve for time δm and recurse
                pop!(distinct_entries)
                pop!(multiplicities)
                multiplicities[end] += km
                distinct_entries[1] -= δm/kp
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


# With only 2 entries, the next crossing isn't just determined by
# gap = min(δp, δm).
function _majmin_onenorm_2entries!(spvm::SortedProbVecMult, ϵ)
    @unpack distinct_entries, multiplicities = spvm
    @assert length(distinct_entries) == length(multiplicities) == 2
    kp, km = multiplicities
    μp, μm = distinct_entries
    # the larger entry is decreasing at a rate t/kp
    # the smaller entry is increasing at a rate t/km
    # Need to solve   `μp - t/kp == μm + t/km` for `t`
    # t = (μp - μm) / ( 1/km + 1/kp )
    t = (μp - μm)*(km*kp) / ( km + kp )
    if ϵ < t
        distinct_entries[1] -= ϵ/kp
        distinct_entries[2] += ϵ/km
        return spvm
    else
        pop!(distinct_entries)
        pop!(multiplicities)
        distinct_entries[] -= t/kp
        multiplicities[] += km
        return spvm
    end
end
