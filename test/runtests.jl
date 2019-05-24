using MajorizationExtrema
using Test, LinearAlgebra

function check_simplexpt(q)
    if eltype(q) <: Rational
        @test sum(q) == one(eltype(q))
    else
        @test sum(q) ≈ one(eltype(q))
    end
    @test all(x -> x >= 0, q)
end

function entropy(q)
    eta(x) = x <= 0 ? zero(x) : -x*log2(x)
    return sum(eta, q)
end

# this likely doesn't sample uniformly
function sample_L1(d)
    u = randsimplexpt(2d)
    return [ u[i] - u[i+d] for i = 1:d ]
end

# this definitely does not sample uniformly
function sample_TV(p, ϵ)
    v = sample_L1(length(p))
    q = p + ϵ*v
    @. q = clamp(q, 0, 1)
    q .= q ./ sum(q)

    # Not sure if this is necessary, but a priori it is
    TV(p, q) > ϵ && return sample_TV(p, ϵ)

    return q
end

@testset "Test local bound" begin
    for d = [2, 3, 5, 20]
        p = randsimplexpt(d)
        entropy_p = entropy(p)
        for ϵ = [.01, .01, .1]
            LB = localbound(entropy, p, ϵ)
            p_max = majmax(p, ϵ)
            p_min = majmin(p, ϵ)
            @test ( LB ≈ entropy(p) - entropy(majmax(p, ϵ)) ) || ( LB ≈ entropy(majmin(p, ϵ)) - entropy(p) )
            for _ = 1:5
                r = sample_TV(p, ϵ)
                @test p_min ≺ r
                @test r ≺ p_max
                @test abs(entropy(r) - entropy_p) <= LB
            end
        end
    end
end

@testset "Basic checks" begin

    p = [1//3, 1//3, 1//3]
    q = [1//4, 1//4, 1//2]

    @test p ≺ q
    @test TV(p, q) == 1//6

    @test majmin(q, 1//20) == [11//40, 11//40, 9//20]
    @test majmax(q, 1//20) == [1//5, 1//4, 11//20]
    @test majmax(q, 2) == [0//1, 0//1, 1//1]
    @test majmin(q, 2) == [1//3, 1//3, 1//3]

    p = float.(p)
    q = float.(q)

    @test p ≺ q
    @test majmin(q, 1/20) ≈ [11/40, 11/40, 9/20]
    @test majmax(q, 1/20) ≈ [1/5, 1/4, 11/20]
    @test majmax(q, 2) ≈ [0, 0, 1]


end

@testset "Tests with type $T" for T = [Float64, Rational{Int}, Rational{BigInt}] 
    for d in [2, 3, 5, 10]
    
        δ = zeros(T, d)
        δ[1] = one(T)

        u = ones(T,d)/T(d)

        for ϵ in [T(1//30), T(1//20), T(1//4)]
            if T == Float64
                r = randsimplexpt(d)
                tol = 1e-8
            elseif T == Rational{Int}
                r = randsimplexpt(d, rand(10:10000))
                tol = zero(T)
            elseif T == Rational{BigInt}
                r = randsimplexpt(d, big(rand(10:10000)))
                tol = zero(T)
            end

            r_max = majmax(r, ϵ)
            r_min = majmin(r, ϵ)

            @test TV(r, r_max) ≈ .5*norm(r-r_max, 1)
            @test TV(r, r_min) ≈ .5*norm(r-r_min, 1)

            @test TV(r, r_max) <= ϵ + tol
            @test TV(r, r_min) <= ϵ + tol

            check_simplexpt(r)
            check_simplexpt(r_max)
            check_simplexpt(r_min)

            @test eltype(r) == T
            @test eltype(r_max) == T
            @test eltype(r_min) == T

            @test r ≺ δ
            @test r ≺ r_max
            @test r_min ≺ r
            @test u ≺ r
        end
    end
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
