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

q = [1//4, 1//4, 1//2]
majmin(q, 1//20) # [11//40, 11//40, 9//20]
majmax(q, 1//20) # [1//5, 1//4, 11//20]
```