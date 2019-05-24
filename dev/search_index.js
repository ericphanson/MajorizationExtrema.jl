var documenterSearchIndex = {"docs":
[{"location":"#MajorizationExtrema.jl-1","page":"Home","title":"MajorizationExtrema.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [MajorizationExtrema]","category":"page"},{"location":"#MajorizationExtrema.:≺-Union{Tuple{T2}, Tuple{T1}, Tuple{AbstractArray{T1,1},AbstractArray{T2,1}}} where T2 where T1","page":"Home","title":"MajorizationExtrema.:≺","text":"≺(p::AbstractVector{T1}, q::AbstractVector{T2}; tol = (T1 <: AbstractFloat || T2 <: AbstractFloat) ? T2(1e-8) : zero(T2) ) where {T1, T2}\n\nReturns true if q majorizes p and false otherwise. The keyword argument tol specifies a tolerance for the comparisons. Can be used in infix form, i.e. p ≺ q.\n\n\n\n\n\n","category":"method"},{"location":"#MajorizationExtrema.TV-Tuple{Any,Any}","page":"Home","title":"MajorizationExtrema.TV","text":"TV(p, q)\n\nReturns the total variation distance between p and q (half the sum of the absolute value of the differences).\n\n\n\n\n\n","category":"method"},{"location":"#MajorizationExtrema.localbound-Tuple{Any,Any,Any}","page":"Home","title":"MajorizationExtrema.localbound","text":"localbound(f, p, ϵ)\n\nReturns δ such that |f(p) - f(q)| <= δ for any q with TV(p, q) <= ϵ, where p and q are probability vectors. Requires f to be Schur convex or Schur concave (but does not check this condition).\n\n\n\n\n\n","category":"method"},{"location":"#MajorizationExtrema.majmax-Tuple{Any,Any}","page":"Home","title":"MajorizationExtrema.majmax","text":"majmax(q, ϵ)\n\nReturns the maximum in majorization order over the total variation ball (of probability vectors) of radius ϵ around a probability vector q.\n\n\n\n\n\n","category":"method"},{"location":"#MajorizationExtrema.majmin-Union{Tuple{ϵT}, Tuple{Any,ϵT}} where ϵT","page":"Home","title":"MajorizationExtrema.majmin","text":"majmin(q, ϵ)\n\nReturns the minimum in majorization order over the total variation ball (of probability vectors) of radius ϵ around a probability vector q.\n\n\n\n\n\n","category":"method"},{"location":"#MajorizationExtrema.randsimplexpt-Tuple{Any}","page":"Home","title":"MajorizationExtrema.randsimplexpt","text":"randsimplexpt(d)\n\nGenerates points uniformly at random on the standard d-1 dimensional simplex using an algorithm by Smith and Tromble.\n\n\n\n\n\n","category":"method"},{"location":"#MajorizationExtrema.randsimplexpt-Union{Tuple{T}, Tuple{Any,T}} where T<:Integer","page":"Home","title":"MajorizationExtrema.randsimplexpt","text":"randsimplexpt(d, N::T) where {T <: Integer}\n\nGenerates rational numbers with denominator at most N uniformly at random on the d-1 dimensional simplex using an algorithm by Smith and Tromble.\n\n\n\n\n\n","category":"method"},{"location":"#MajorizationExtrema.simplexpt-Tuple{Any}","page":"Home","title":"MajorizationExtrema.simplexpt","text":"simplexpt(unif)\n\nTakes a vector of length d-1 of numbers between 0 and 1 and converts it a point on the standard d dimensional simplex.\n\n\n\n\n\n","category":"method"}]
}