include("fisher_yates!.jl")

"""
    simplexpt(unif)

Takes a vector of length `d-1` of numbers between `0` and `1` and converts it a point on the standard `d-1` dimensional simplex.
"""
function simplexpt(unif)
    d = length(unif) + 1; T = eltype(unif)
    w= zeros(T, d+1)
    w[2:d] .= sort(unif)
    w[d+1] = one(T)
    diff(w)
end

"""
    randprobvec([rng], d::T) where {T <: Integer}

Generates a random probability vector of length `d` using an algorithm by [Smith and Tromble](http://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf).
"""
randprobvec(rng::AbstractRNG, d::T) where {T <: Integer} = simplexpt(rand(rng, d-1))

randprobvec(d::T) where {T <: Integer} = randprobvec(Random.GLOBAL_RNG, d)

"""
    randprobvec([rng], d::dT, N::T) where {T <: Integer, dT <: Integer}

Generates a random probability vector of length `d` whose entries are rational numbers with denominator at most `N` using an algorithm by [Smith and Tromble](http://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf).
"""
function randprobvec(rng::AbstractRNG, d::dT, N::T) where {T <: Integer, dT <: Integer}
    unif = zeros(T, d-1)
    fisher_yates_sample!(rng, 1:N, unif)
    unif = unif .// N
    simplexpt(unif)
end

randprobvec(d::dT, N::T) where {T <: Integer, dT <: Integer} = randprobvec(Random.GLOBAL_RNG, d, N)


"""
    randunitary([rng], d::T) where {T <: Integer}

Generates a unitary matrix of dimension `d` at random according to the Haar measure, using an algorithm described by Maris Ozols in ["How to generate a random unitary matrix"](http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20%5Bpaper%5D.pdf).
"""
function randunitary(rng::AbstractRNG, d::T) where {T <: Integer}
   RG = randn(rng, d, d) + im*randn(rng, d, d)
   Q,R = qr!(RG);
   r = diag(R)
   L = Diagonal(r./abs.(r));
   return Q*L
end

randunitary(d::T) where {T <: Integer} = randunitary(Random.GLOBAL_RNG, d)


"""
    randdm([rng], d::T) where {T <: Integer}

Generates a density matrix of dimension `d` at random.
"""
function randdm(rng, d::T) where {T <: Integer}
    eigs = Diagonal(randprobvec(rng, d))
    U = randunitary(rng, d)
    ρ =  U * eigs * (U')
    return Hermitian((ρ + ρ')/2)
end

randdm(d::T) where {T <: Integer} = randdm(Random.GLOBAL_RNG, d)
