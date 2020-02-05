module MajorizationExtrema

using LinearAlgebra, Random
using UnPack
using LinearAlgebra: dot
using ForwardDiff

export TV, tracedist, ≺, ≻ # distances, majorization
export majmax, majmin, localbound # majorization-extrema
export randprobvec, randunitary, randdm # generate random data


include("random_utilities.jl")

include("classical_case.jl")

include("probvecmults.jl")

# reduce to classical via eigenvalues
include("quantum_case.jl")

include("entropies.jl")


"""
    localbound(f, p, ϵ)

Requires `f` to be Schur convex or Schur concave (but does not check this condition). 

* If `p` is a vector, it is a assumed to be a probability vector, and `localbound` returns `δ` such that `|f(p) - f(q)| <= δ` for any probability vector `q` with `TV(p, q) <= ϵ`.
* If `p` is a matrix, it is assumed to be a density matrix, and `localbound` returns `δ` such that `|f(p) - f(q)| <= δ` for any density matrix `q` with `tracedist(p, q) <= ϵ`.
"""
function localbound(f, p, ϵ)
    fp = f(p)
    f_max = f(majmax(p, ϵ))
    f_min = f(majmin(p, ϵ))
    return max(abs(fp - f_max), abs(fp - f_min))
end

end # module
