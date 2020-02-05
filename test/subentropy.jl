using MajorizationExtrema: subentropy, entropy

function _naive_subentropy(r::AbstractVector)
    @assert length(unique(r)) == length(r)
    n = length(r)
    f = x -> x^n * log(x)
    s = zero(eltype(r))
    for i in eachindex(r)
        d = prod(r[i] - r[j] for j = eachindex(r) if j != i)
        s -= r[i]^n * log(r[i]) / d
    end
    return s
end

@testset "Subentropy" begin
    for T in (Float64, BigFloat), d in (2, 3, 5)
        r = randprobvec(T, d)
        length(unique(r)) == d || continue
        @test _naive_subentropy(r) ≈ subentropy(r)

        # basic bounds
        @test 0 <= subentropy(r) <= entropy(r)
    end
end
