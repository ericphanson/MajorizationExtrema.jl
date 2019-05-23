using MajorizationExtrema
using Test

function check_simplexpt(q)
    if eltype(q) <: Rational
        @test sum(q) == one(eltype(q))
    else
        @test sum(q) ≈ one(eltype(q))
    end
    @test all(x -> x >= 0, q)
end

@testset "Basic checks" begin

    p = [1//3,1//3,1//3]
    q = [1//4,1//4,1//2]

    @test p ≺ q
    @test TV(p, q) == 1//6

    @test majmin(q, 1//20) == [11//40, 11//40, 9//20]

    @test majmax(q, 1//20) == [1//5, 1//4, 11//20]

    @test majmax(q, 2) == [0//1, 0//1, 1//1]

    @test majmin(q, 2) == [1//3, 1//3, 1//3]


end

@testset "Random tests" begin
    for _ = 1:10
        q_rational = randsimplexpt(3, 1000)
        check_simplexpt(q_rational)

        q_float = float.(q_rational)
        check_simplexpt(q_float)

        for ϵ in (1//100, 2//100)
            @test majmin(q_rational, ϵ) ≈ majmin(q_float, ϵ)
            @test majmax(q_rational, ϵ) ≈ majmax(q_float, ϵ)
        end

        for ϵ in (1, 2, 2//3)
            @test majmin(q_rational, ϵ) == [1//3, 1//3, 1//3]
            @test majmin(q_float, ϵ) ≈ [1/3, 1/3, 1/3]

            @test sort(majmax(q_rational, ϵ)) == [0//1, 0//1, 1//1]
            @test sort(majmax(q_float, ϵ)) ≈ [0.0, 0.0, 1.0]
        end
    end
end
