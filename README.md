# MajorizationExtrema

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ericphanson.github.io/MajorizationExtrema.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ericphanson.github.io/MajorizationExtrema.jl/dev)
[![Build Status](https://travis-ci.com/ericphanson/MajorizationExtrema.jl.svg?branch=master)](https://travis-ci.com/ericphanson/MajorizationExtrema.jl)
[![Codecov](https://codecov.io/gh/ericphanson/MajorizationExtrema.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ericphanson/MajorizationExtrema.jl)

## Install

This package isn't registered yet, but you can install it via

```julia
] add https://github.com/ericphanson/MajorizationExtrema.jl
```

## Usage

Compute the majorization minimizer and maximizer over total variation balls (in finite dimensions).

```julia
using MajorizationExtrema

q = [0.25, 0.25, 0.5]
majmin(q, 0.05) # [0.275, 0.275, 0.45]
majmax(q, 0.05) # [0.2, 0.25, 0.55]
```

Also works with rational numbers:

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

Based on:

>E. P. Hanson and N. Datta, “Maximum and minimum entropy states yielding local continuity bounds,” Journal of Mathematical Physics, vol. 59, no. 4, p. 042204, Apr. 2018.

for which a preprint is available here: <https://arxiv.org/abs/1706.02212v2>.