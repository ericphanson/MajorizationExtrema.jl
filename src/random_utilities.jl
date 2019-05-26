include("fisher_yates!.jl")

"""
    simplexpt(unif)

Takes a vector of length `d-1` of numbers between `0` and `1` and converts it a point on the standard `d` dimensional simplex.
"""
function simplexpt(unif)
    d = length(unif) + 1; T = eltype(unif)
    w= zeros(T,d+1)
    w[2:d] .= sort(unif)
    w[d+1] = one(T)
    diff(w)
end

"""
    randsimplexpt(d)

Generates points uniformly at random on the standard `d-1` dimensional simplex using an algorithm by [Smith and Tromble](http://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf).
"""
randsimplexpt(d) = simplexpt(rand(d-1))

"""
    randsimplexpt(d, N::T) where {T <: Integer}

Generates rational numbers with denominator at most `N` uniformly at random on the `d-1` dimensional simplex using an algorithm by [Smith and Tromble](http://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf).
"""
function randsimplexpt(d, N::T) where {T <: Integer}
    unif = zeros(T, d-1)
    fisher_yates_sample!(1:N, unif)
    unif = unif .// N
    simplexpt(unif)
end


"""
    randunitary(d)

Generates a unitary matrix of dimension `d` at random according to the Haar measure, using an algorithm described by Maris Ozols in ["How to generate a random unitary matrix"](http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20%5Bpaper%5D.pdf).
"""
function randunitary(d)
   RG = randn(d,d) + im*randn(d,d)
   Q,R = qr!(RG);
   r = diag(R)
   L = Diagonal(r./abs.(r));
   return Q*L
end

"""
    randdm(d)

Generates a density matrix of dimension `d` at random.
"""
function randdm(d)
    eigs = Diagonal(randsimplexpt(d))
    U = randunitary(d)
    ρ =  U * eigs * (U')
    return Hermitian((ρ + ρ')/2)
end