# MajorizationExtrema

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ericphanson.github.io/MajorizationExtrema.jl/dev)
[![Build Status](https://travis-ci.com/ericphanson/MajorizationExtrema.jl.svg?branch=master)](https://travis-ci.com/ericphanson/MajorizationExtrema.jl)
[![Codecov](https://codecov.io/gh/ericphanson/MajorizationExtrema.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ericphanson/MajorizationExtrema.jl)

## Install

This package isn't registered yet, but you can install it via

```julia
] add https://github.com/ericphanson/MajorizationExtrema.jl
```

This package implements methods from

>E. P. Hanson and N. Datta, “Maximum and minimum entropy states yielding local continuity bounds,” Journal of Mathematical Physics, vol. 59, no. 4, p. 042204, Apr. 2018.

for which a preprint is available here: <https://arxiv.org/abs/1706.02212v2>.

## Usage

Provides the following functions:

```julia
TV, tracedist, ≺ # distances, majorization
majmax, majmin, localbound # majorization-extrema
randprobvec, randunitary, randdm # generate random data
```

Compute the majorization minimizer and maximizer over total variation balls of discrete probability distributions (in finite dimensions).

```julia
using MajorizationExtrema

q = [0.25, 0.25, 0.5]
majmin(q, 0.05) # [0.275, 0.275, 0.45]
majmax(q, 0.05) # [0.2, 0.25, 0.55]
```

These vectors satisfy `majmin(q, 0.05) ≺ p ≺ majmax(q, 0.05)` for any probability distribution `p` within `0.05` of `q` in total variation distance (`TV(p, q) <= 0.05`), where `≺` is the [majorization preorder](https://en.wikipedia.org/wiki/Majorization).

These functions also work with rational numbers:

```julia
q = [1//4, 1//4, 1//2]
majmin(q, 1//20) # [11//40, 11//40, 9//20]
majmax(q, 1//20) # [1//5, 1//4, 11//20]
```

If you're doing extended computations with rational numbers, you may want to use `Rational{BigInt}` types which won't overflow:

```julia
q = big.([1//4, 1//4, 1//2])
...
```

The existence of majorization minimizer and maximizers allows the computation of local bounds of Schur concave or Schur convex functions on probability distributions around a specific distribution; `localbound` is provided as a convenience function for this. For example,

```julia
using MajorizationExtrema

function entropy(q)
    eta(x) = x <= 0 ? zero(x) : -x*log2(x)
    return sum(eta, q)
end

p = [.2, .3, .5] # a fixed probability distribution
q = randprobvec(3) # a random probability distribution
ϵ = TV(p, q) # total variation distance
abs(entropy(p) - entropy(q)) <= localbound(entropy, p, ϵ) # true
```

Note `localbound(entropy, p, ϵ)` only depends on the choice of `q` via `ϵ`; the bound holds uniformly over all probability distributions `q` within total variation distance `ϵ` of `p`.

These functions also work on quantum density matrices. For example,

```julia
julia> using MajorizationExtrema

julia> ρ = randdm(3)
3×3 LinearAlgebra.Hermitian{Complex{Float64},Array{Complex{Float64},2}}:
   0.235005+0.0im         0.00492117+0.00548039im  -0.0168253-0.020445im 
 0.00492117-0.00548039im    0.458451+0.0im         -0.0724026+0.0714849im
 -0.0168253+0.020445im    -0.0724026-0.0714849im     0.306544+0.0im      

julia> σ = majmin(ρ, .01)
3×3 LinearAlgebra.Hermitian{Complex{Float64},Array{Complex{Float64},2}}:
   0.242182+0.0im         0.00433537+0.00652348im  -0.0142623-0.0168348im
 0.00433537-0.00652348im     0.45091+0.0im         -0.0689234+0.0679345im
 -0.0142623+0.0168348im   -0.0689234-0.0679345im     0.306909+0.0im      

julia> σ ≺ ρ
true
```

(One will likely want to write `using LinearAlgebra` when working with quantum states in order to compute eigenvalues, etc.)