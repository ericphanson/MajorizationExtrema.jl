using CairoMakie, MajorizationExtrema
using MajorizationExtrema: InfNorm, SortedProbVecMult
using LinearAlgebra


# this likely doesn't sample uniformly
function sample_L1(d)
    u = randprobvec(2d)
    return [u[i] - u[i+d] for i = 1:d]
end

# this definitely does not sample uniformly
function sample_TV(p, ϵ)
    v = sample_L1(length(p))
    q = @. clamp(p + ϵ * v, 0, 1)
    q .= q ./ sum(q)

    # Not sure if this is necessary, but a priori it is
    TV(p, q) > ϵ && return sample_TV(p, ϵ)

    return q
end

function sample_inf(d)
    v = 2rand(d) .- 1
    return v .- sum(v)
end

function sample_inf(p, ϵ)
   v = sample_inf(length(p))
   q = @. clamp(p + ϵ * v, zero(eltype(p)), one(eltype(p)))
   q .= q ./ sum(q)
   any(isnan, q) && return sample_inf(p, ϵ)
   norm(p - q, Inf) > ϵ && return sample_inf(p, ϵ)

   q = @. rationalize(BigInt, q)
   @. q = clamp(q, zero(eltype(q)), one(eltype(q)))
   q .= q .// sum(q)
   norm(p - q, Inf) > ϵ && return sample_inf(p, ϵ)
   return q
end


function lorenz!(ax, p; kw...)
    p = sort(p; rev=false)
    x = collect((0:length(p)) ./ length(p))
    y = [0; cumsum(p)] ./ sum(p)
    x = Float64.(x)
    y = Float64.(y)
    scatterlines!(ax, x, y; kw...)
end


# Makie.convert_arguments(::Type{<:Scatter}, L::Lorenz) = (collect((1:length(L.p)) ./ length(L.p)), cumsum(sort(L.p; rev=true)) ./ sum(L.p))
# plottype(::Lorenz) = ScatterLines

d = 5

function example(d, ϵ; N = 100, p)
    q = randprobvec(d, big(100_000))
    ϵ = big(rationalize(ϵ))
    sort!(q; rev=true)
    if p == 1
        nrm = OneNorm()
        sample = () -> sample_TV(q, ϵ)
    elseif p == Inf
        nrm = InfNorm()
        sample = () -> sample_inf(q, ϵ)
    else
        error("`p` must be 1 or Inf; got $p")
    end
    fig = Figure()
    ax = Axis(fig[1,1])
    lorenz!(ax, q; label="q")

    q_min = majmin(nrm, SortedProbVecMult(q), ϵ)

    lorenz!(ax,q_min; label="min")

    for _ = 1:N
        r = sample()
        @assert norm(q - r, p) <= ϵ norm(q - r, p)
        @assert sum(r) ≈ 1
        @assert all(>=(0), r)
        lorenz!(ax, r, color=(:gray, 0.65), markersize=2, linewidth=0.5)
        if !(q_min ≺ r)
            return (; q_min, r, q)
        end
    end
    axislegend(ax; position=:lt)
    fig
end


# Counterexample

r = BigFloat[0.5029277026173153953901596916818025377675264796549552335054911802433572658566893, 0.2532618958883083476845783210120004920380552411564360036445550734061488887917368, 0.09730765240217434870629777804185540217187642157964938793634674296073053856139903, 0.06755518200020205877091445188399030380035851678647313594812135305483526594141762, 0.07894756709199984944804975738035126422218334082248623896548565033492804084876043]
q_min = Rational{BigInt}[12447//25000, 737//2500, 5183//75000, 5183//75000, 5183//75000]
q = Rational{BigInt}[14947//25000, 17117//50000, 4693//100000, 877//100000, 51//12500]

sort!(r; rev=true)
sort!(q; rev=true)
sort!(q_min; rev=true)
r= @. rationalize(BigInt, r)
@. r = clamp(r, zero(eltype(r)), one(eltype(r)))
r .= r .// sum(r)

norm(r - q, Inf) < 0.1 # true
q_min ≺ q # true
q_min ≺ r # false


q_min_real = MajorizationExtrema.majmin_inf(q, 0.1)
# julia> float(q_min)
# 5-element Vector{BigFloat}:
#  0.4978800000000000000000000000000000000000000000000000000000000000000000000000014
#  0.2947999999999999999999999999999999999999999999999999999999999999999999999999989
#  0.06910666666666666666666666666666666666666666666666666666666666666666666666666658
#  0.06910666666666666666666666666666666666666666666666666666666666666666666666666658
#  0.06910666666666666666666666666666666666666666666666666666666666666666666666666658

# julia> float(q_min_real)
# 5-element Vector{BigFloat}:
#  0.4978799999999999944488848768742172978815968973143114714876453076884424316663615
#  0.2423399999999999944488848768742172978815857647652187539259532246226679932894977
#  0.08659333333333333703407674875052180141362866746072654555922350826047289788845693
#  0.08659333333333333703407674875052180141206903045483947723165730617671990014629295
#  0.08659333333333333703407674875052180141107400917380593636554503517153825150738743

# Simpler:
δ = 1//100
q = [1//2 - δ, 1//2 - δ, δ, δ]
# julia> MajorizationExtrema.majmin_inf(q, 1//10)
# 4-element Vector{BigFloat}:
#  0.3900000000000000133226762955018784850835800169813509748601599783286893520877027
#  0.3900000000000000133226762955018784850835800170867526693544912843507680059225947
#  0.109999999999999986677323704498121514916363282037248558687358814688779246729992
#  0.109999999999999986677323704498121514916422507645402385259588037867006404128714

# julia> float.(majmin(InfNorm(), SortedProbVecMult(q), 1//big(10)))
# (Float64.(distinct_entries), Float64(ϵ), multiplicities) = ([0.49, 0.01], 0.1, [2, 2])
# 4-element Vector{BigFloat}:
#  0.4400000000000000000000000000000000000000000000000000000000000000000000000000014
#  0.4400000000000000000000000000000000000000000000000000000000000000000000000000014
#  0.06000000000000000000000000000000000000000000000000000000000000000000000000000024
#  0.06000000000000000000000000000000000000000000000000000000000000000000000000000024

# Notes:
# `MajorizationExtrema.majmin_inf` actually computes the majorization infimum
# but that infimum seems to be in the ball, which makes it actually a min.
