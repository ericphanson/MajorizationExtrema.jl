include("fisher_yates!.jl")

"""
    simplexpt(unif)

Takes a vector of length `d-1` of numbers between `0` and `1` and converts it a point on the standard `d-1` dimensional simplex.
"""
function simplexpt(unif)
    d = length(unif) + 1; T = eltype(unif)
    w = zeros(T, d + 1)
    w[2:d] .= sort(unif)
    w[d + 1] = one(T)
    diff(w)
end

"""
    randprobvec([rng], d::Int)

Generates a random probability vector of length `d` using an algorithm by [Smith and Tromble](http://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf).
"""
randprobvec(rng::AbstractRNG, ::Type{N},  d::Int) where {N <: Number} = simplexpt(rand(rng, N, d - 1))

randprobvec(rng::AbstractRNG, d::Int) = simplexpt(rand(rng, Float64, d - 1))

randprobvec(::Type{N}, d::Int) where {N <: Number} = randprobvec(Random.GLOBAL_RNG, N, d)

randprobvec(d::Int) = randprobvec(Random.GLOBAL_RNG, Float64, d)


"""
    randprobvec([rng], d::Int, N::T) where {T <: Integer}

Generates a random probability vector of length `d` whose entries are rational numbers with denominator at most `N` using an algorithm by [Smith and Tromble](http://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf).
"""
function randprobvec(rng::AbstractRNG, d::Int, N::T) where {T <: Integer}
    unif = zeros(T, d - 1)
    fisher_yates_sample!(rng, 1:N, unif)
    unif = unif .// N
    simplexpt(unif)
end

randprobvec(d::Int, N::T) where {T <: Integer} = randprobvec(Random.GLOBAL_RNG, d, N)


"""
    randunitary([rng], d::Int)

Generates a unitary matrix of dimension `d` at random according to the Haar measure, using an algorithm described by Maris Ozols in ["How to generate a random unitary matrix"](http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20%5Bpaper%5D.pdf).
"""
function randunitary(rng::AbstractRNG, d::Int)
    RG = randn(rng, d, d) + im * randn(rng, d, d)
    Q, R = qr!(RG);
    r = diag(R)
    L = Diagonal(r ./ abs.(r));
    return Q * L
end

randunitary(d::Int)  = randunitary(Random.GLOBAL_RNG, d)


"""
    randdm([rng], d::Int)

Generates a density matrix of dimension `d` at random.
"""
function randdm(rng, d::Int)
    eigs = Diagonal(randprobvec(rng, d))
    U = randunitary(rng, d)
    ρ =  U * eigs * (U')
    return Hermitian((ρ + ρ') / 2)
end

randdm(d::Int) = randdm(Random.GLOBAL_RNG, d)

"""
    randmajpair([rng], [T], d::Int, [N])

Generates a pair of probability vectors `p`, `q` of length `d`
such that `p ≺ q` by rejection sampling. Supply a type `T` to
choose the numeric type of the elements of the vectors, and an
integer `N` to use rational numbers with maximum denominator `N`.
"""
function randmajpair(args...)
    while true
        p = randprobvec(args...)
        q = randprobvec(args...)
        p ≺ q && return p, q
        p ≻ q && return q, p
    end
end

"""
    randincomppair([rng], [T], d::Int, [N]; sortby = entropy)

Generates a pair of probability vectors `p`, `q` of length `d`
such that `!(p ≺ q)` and `!(q ≺ p)` by rejection sampling. Supply
a type `T` to choose the numeric type of the elements of the vectors,
and an integer `N` to use rational numbers with maximum denominator `N`.

Optionally, pass a function (from probability vectors to an
ordered set) to the keyword argument `sortby` to sort the results.
"""
function randincomppair(args...; sortby = entropy)
    while true
        p = randprobvec(args...)
        q = randprobvec(args...)
        p ≺ q && continue
        p ≻ q && continue
        if entropy(p) >= entropy(q)
            return p, q
        else
            return q, p
        end
    end
end