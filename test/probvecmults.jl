@testset "ProbVecMult roundtrip" begin
    for d in (2, 3, 5), T in (Float64, BigFloat)
        v = randprobvec(T, d)
        v2 = copy(v)
        result = collect(SortedProbVecMult(v))
        @test v2 == v # not modified
        @test issorted(result; rev = true) # sorted
        @test result ≈ sort(v2; rev = true) # roundtrip successfully
    end
end

@testset "ProbVecMult `majmin`" begin

    for d in (2, 3, 5), T in (Float64, BigFloat)
        v = randprobvec(T, d)
        spvm = SortedProbVecMult(v)
        for ϵ in range(0.0, 1.0, length = 100)
            @test majmin(sort(v; rev = true), ϵ) ≈ collect(majmin(spvm, ϵ))
        end
    end

end
