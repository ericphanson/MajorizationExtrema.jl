# The following is taken from StatsBase
# https://github.com/JuliaStats/StatsBase.jl/blob/master/src/sampling.jl
# which is available under the following MIT license:

# Copyright (c) 2012-2016: Dahua Lin, Simon Byrne, Andreas Noack, Douglas Bates, John Myles White, Simon Kornblith, and other contributors.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


"""
    fisher_yates_sample!([rng], a::AbstractArray, x::AbstractArray)
Fisher-Yates shuffling (with early termination).
Pseudo-code:
```
n = length(a)
k = length(x)
# Create an array of the indices
inds = collect(1:n)
for i = 1:k
    # swap element `i` with another random element in inds[i:n]
    # set element `i` in `x`
end
```
This algorithm consumes `k=length(x)` random numbers. It uses an integer array of
length `n=length(a)` internally to maintain the shuffled indices. It is considerably
faster than Knuth's algorithm especially when `n` is greater than `k`.
It is ``O(n)`` for initialization, plus ``O(k)`` for random shuffling
"""
function fisher_yates_sample!(rng::AbstractRNG, a::AbstractArray, x::AbstractArray)
    n = length(a)
    k = length(x)
    k <= n || error("length(x) should not exceed length(a)")

    inds = Vector{Int}(undef, n)
    for i = 1:n
        @inbounds inds[i] = i
    end

    @inbounds for i = 1:k
        j = rand(rng, i:n)
        t = inds[j]
        inds[j] = inds[i]
        inds[i] = t
        x[i] = a[t]
    end
    return x
end
fisher_yates_sample!(a::AbstractArray, x::AbstractArray) =
fisher_yates_sample!(Random.GLOBAL_RNG, a, x)