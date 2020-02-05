using MajorizationExtrema
using MajorizationExtrema: entropy
using Test, LinearAlgebra

function check_simplexpt(q)
    if eltype(q) <: Rational
        @test sum(q) == one(eltype(q))
    else
        @test sum(q) ≈ one(eltype(q))
    end
    @test all(x -> x >= 0, q)
end

function check_dm(ρ)
    @test tr(ρ) ≈ one(eltype(ρ))
    @test eigmin(ρ) >= -1e-8
end

# this likely doesn't sample uniformly
function sample_L1(d)
    u = randprobvec(2d)
    return [u[i] - u[i+d] for i = 1:d]
end

# this definitely does not sample uniformly
function sample_TV(p, ϵ)
    v = sample_L1(length(p))
    q = p + ϵ * v
    @. q = clamp(q, 0, 1)
    q .= q ./ sum(q)

    # Not sure if this is necessary, but a priori it is
    TV(p, q) > ϵ && return sample_TV(p, ϵ)

    return q
end

@testset "Types" begin
    @test @inferred(randdm(3)) isa Hermitian
    ρ = randdm(3)
    @test @inferred(majmin(ρ, 0.01)) isa Hermitian
    @test @inferred(majmax(ρ, 0.01)) isa Hermitian
    @test eltype(@inferred(majmin(ρ, 1 // 10))) == Complex{Float64}

    @test typeof(@inferred(randprobvec(3))) == Vector{Float64}
    @test typeof(@inferred(randprobvec(3, 100))) == Vector{Rational{Int}}
    @test typeof(@inferred(randprobvec(3, big(100)))) == Vector{Rational{BigInt}}
    @test typeof(@inferred(randprobvec(BigFloat, 3))) == Vector{BigFloat}

    for m in (majmin, majmax)
        p = randprobvec(3)
        @test typeof(@inferred(m(p, 0.01))) == Vector{Float64}
        @test typeof(@inferred(m(p, 1 // 10))) == Vector{Float64}
        @test typeof(@inferred(m(p, big(1 // 10)))) == Vector{BigFloat}

        p = randprobvec(3, 100)
        @test typeof(@inferred(m(p, 0.01))) == Vector{Float64}
        @test typeof(@inferred(m(p, 1 // 10))) == Vector{Rational{Int}}
        @test typeof(@inferred(m(p, big(1 // 10)))) == Vector{Rational{BigInt}}

        p = randprobvec(3, big(100))
        @test typeof(@inferred(m(p, 0.01))) == Vector{BigFloat}
        @test typeof(@inferred(m(p, 1 // 10))) == Vector{Rational{BigInt}}
        @test typeof(@inferred(m(p, big(1 // 10)))) == Vector{Rational{BigInt}}
    end

end

@testset "Quantum states" begin

    X = [1.0 2.0; 0.0 0.0]
    @test_throws ArgumentError majmin(X, 0.01)
    @test_throws ArgumentError majmax(X, 0.01)

    ρ = randdm(2)
    @test_throws ArgumentError ρ ≺ X
    @test_throws ArgumentError X ≺ ρ
    for d in (2, 3, 5)
        idmat = Matrix(1.0I, d, d)

        for ϵ in (1, 2, 1 - 1 / d + 0.01)
            ρ = randdm(d)
            check_dm(ρ)
            @test majmin(ρ, ϵ) ≈ idmat / d
            @test tr(majmax(ρ, ϵ)^2) ≈ 1.0
            @test majmin(ρ, ϵ) ≺ ρ
            @test ρ ≺ majmax(ρ, ϵ)
        end

        for ϵ in (0.01, 0.05, 0.1)
            ρ = randdm(d)
            check_dm(ρ)
            @test localbound(entropy, ρ, ϵ) <= ϵ * log2(d - 1) + entropy([ϵ, 1 - ϵ])

            p = randprobvec(d)
            U = randunitary(d)
            @test U * U' ≈ idmat
            @test U' * U ≈ idmat

            ρ = U * Diagonal(p) * U'
            ρ = (ρ + ρ') / 2
            check_dm(ρ)

            @test majmin(ρ, ϵ) ≈ U * Diagonal(majmin(p, ϵ)) * U'
            @test majmax(ρ, ϵ) ≈ U * Diagonal(majmax(p, ϵ)) * U'

            @test majmin(ρ, ϵ) ≺ ρ
            @test ρ ≺ majmax(ρ, ϵ)

            if tracedist(ρ, idmat / d) > ϵ
                @test tracedist(majmin(ρ, ϵ), ρ) ≈ ϵ
            else
                @test majmin(ρ, ϵ) ≈ idmat / d
            end

            p = randprobvec(d)
            q = randprobvec(d)
            ρ = U * Diagonal(p) * U'
            σ = U * Diagonal(q) * U'
            @test tracedist(ρ, σ) ≈ tracedist(Hermitian(ρ), Hermitian(σ))
            @test tracedist(ρ, σ) ≈ TV(p, q)
        end

    end
end

@testset "Test local bound" begin
    for d in [2, 3, 5, 20]
        p = randprobvec(d)
        entropy_p = entropy(p)
        for ϵ in [0.01, 0.01, 0.1]
            LB = localbound(entropy, p, ϵ)
            p_max = majmax(p, ϵ)
            p_min = majmin(p, ϵ)
            @test (LB ≈ entropy(p) - entropy(majmax(p, ϵ))) ||
                  (LB ≈ entropy(majmin(p, ϵ)) - entropy(p))
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

    p1 = [1 // 3, 1 // 3, 1 // 3]
    p2 = [1 // 3, 1 // 3, 1 // 3, 1 // 3]
    @test_throws DimensionMismatch p1 ≺ p2
    @test_throws DimensionMismatch float.(p1) ≺ float.(p2)
    @test_throws DimensionMismatch TV(p1, p2)
    @test_throws DimensionMismatch TV(float.(p1), float.(p2))

    p = [1 // 3, 1 // 3, 1 // 3]
    q = [1 // 4, 1 // 4, 1 // 2]

    @test p ≺ q
    @test !(q ≺ p)
    @test TV(p, q) == 1 // 6

    @test majmin(q, 1 // 20) == [11 // 40, 11 // 40, 9 // 20]
    @test majmax(q, 1 // 20) == [1 // 5, 1 // 4, 11 // 20]
    @test majmax(q, 2) == [0 // 1, 0 // 1, 1 // 1]
    @test majmin(q, 2) == [1 // 3, 1 // 3, 1 // 3]

    p = float.(p)
    q = float.(q)

    @test p ≺ q
    @test !(q ≺ p)
    @test majmin(q, 1 / 20) ≈ [11 / 40, 11 / 40, 9 / 20]
    @test majmax(q, 1 / 20) ≈ [1 / 5, 1 / 4, 11 / 20]
    @test majmax(q, 2) ≈ [0, 0, 1]
    @test majmin(q, 2) ≈ [1 / 3, 1 / 3, 1 / 3]

    @test TV(p, q) ≈ 1 / 6


end

@testset "Tests with type $T" for T in [Float64, BigFloat, Rational{Int}, Rational{BigInt}]
    for d in [2, 3, 5, 10]
        δ = zeros(T, d)
        δ[1] = one(T)

        u = ones(T, d) / T(d)

        for ϵ in [T(1 // 30), T(1 // 20), T(1 // 4)]
            if T == Float64
                r = randprobvec(d)
                tol = 1e-8
            elseif T == BigFloat
                r = randprobvec(BigFloat, d)
                tol = 1e-25
            elseif T == Rational{Int}
                r = randprobvec(d, rand(10:10000))
                tol = zero(T)
            elseif T == Rational{BigInt}
                r = randprobvec(d, big(rand(10:10000)))
                tol = zero(T)
            end

            r_max = majmax(r, ϵ)
            r_min = majmin(r, ϵ)

            @test TV(r, r_max) ≈ 0.5 * norm(r - r_max, 1)
            @test TV(r, r_min) ≈ 0.5 * norm(r - r_min, 1)

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
        q_rational = randprobvec(3, 1000)
        check_simplexpt(q_rational)

        q_float = float.(q_rational)
        check_simplexpt(q_float)

        for ϵ in (1 // 100, 2 // 100)
            @test majmin(q_rational, ϵ) ≈ majmin(q_float, ϵ)
            @test majmax(q_rational, ϵ) ≈ majmax(q_float, ϵ)
        end

        for ϵ in (1, 2, 2 // 3)
            @test majmin(q_rational, ϵ) == [1 // 3, 1 // 3, 1 // 3]
            @test majmin(q_float, ϵ) ≈ [1 / 3, 1 / 3, 1 / 3]

            @test sort(majmax(q_rational, ϵ)) == [0 // 1, 0 // 1, 1 // 1]
            @test sort(majmax(q_float, ϵ)) ≈ [0.0, 0.0, 1.0]
        end

        for T in (Float64, BigFloat)
            p, q = @inferred(MajorizationExtrema.randincomppair(T, 4))
            @test !(p ≺ q) && !(q ≺ p)
            @test !(p ≻ q) && !(q ≻ p)
            check_simplexpt(p)
            check_simplexpt(q)
            @test eltype(p) == T
            @test eltype(q) == T

            p, q = @inferred(MajorizationExtrema.randmajpair(T, 4))
            @test (p ≺ q) && (q ≻ p)
            check_simplexpt(p)
            check_simplexpt(q)
            @test eltype(p) == T
            @test eltype(q) == T
        end


    end
end

include("probvecmults.jl")
include("subentropy.jl")
