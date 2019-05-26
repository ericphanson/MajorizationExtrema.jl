@inline function eig_promote(f, X::Hermitian)
    evals, U = eigen(X)
    return Hermitian(U * Diagonal(f(evals)) * U')
end

"""
    majmin(ρ::Hermitian, ϵ)

Returns the minimum in majorization order over the trace distance ball (of density matrices) of radius `ϵ` around a density matrix `ρ`.
"""
majmin(ρ::Hermitian, ϵ) = eig_promote(p -> majmin(p, ϵ), ρ)

"""
    majmax(ρ::Hermitian, ϵ)

Returns the maximum in majorization order over the trace distance ball (of density matrices) of radius `ϵ` around a density matrix `ρ`.
"""
majmax(ρ::Hermitian, ϵ) = eig_promote(p -> majmax(p, ϵ), ρ)


"""
    ≺(ρ::Hermitian, σ::Hermitian; tol=1e-8)

Returns true if `σ` majorizes `ρ` and false otherwise. The keyword argument `tol` specifies a tolerance for the comparisons. Can be used in infix form, i.e. `ρ ≺ σ`.
"""
≺(ρ::Hermitian, σ::Hermitian; tol=1e-8) = ≺(eigvals(ρ), eigvals(σ); tol=tol)  

"""
    tracedist(ρ::AbstractMatrix, σ::AbstractMatrix)

Computes the trace distance between two matrices `ρ` and `σ`.
"""
tracedist(ρ::AbstractMatrix, σ::AbstractMatrix) = .5*sum(abs, eigvals(ρ - σ))

tracedist(ρ::Hermitian, σ::Hermitian) = .5*sum(svdvals(ρ - σ))

# error checking and wrapping for non `Hermitian`-type matrices:

function ≺(ρ::AbstractMatrix, σ::AbstractMatrix; tol=1e-8)
    ishermitian(ρ) && ishermitian(σ) || throw(ArgumentError("Inputs must be Hermitian matricies"))
    ≺(Hermitian(ρ), Hermitian(σ); tol=tol)  
end

function majmax(ρ::AbstractMatrix, ϵ)
    ishermitian(ρ) || throw(ArgumentError("The input must be a Hermitian matrix"))
    return majmax(Hermitian(ρ), ϵ)
end

function majmin(ρ::AbstractMatrix, ϵ)
    ishermitian(ρ) || throw(ArgumentError("The input must be a Hermitian matrix"))
    return majmin(Hermitian(ρ), ϵ)
end