import LinearAlgebra
using LinearAlgebra: Diagonal
import Random

"""
    random_unitary(n::Integer)

Return a Haar-random nxn unitary matrix.

# Reference
F Mezzadri, arXiv:math-ph/0609050
"""
function random_unitary(n::Integer)
    return n == 2 ? _random_unitary_2x2() : _random_unitary(n)
end

# FIXME: 2x2 allocates a new matrix
function random_unitary!(A::Matrix{ComplexF64})
    dims = size(A)
    if dims[1] != dims[2]
        error("Expecting a square matrix")
    end
    n = dims[1]
    return n == 2 ? _random_unitary_2x2!(A) : _random_unitary!(A)
end

"""
    _random_unitary_2x2()

Return a Haar-random 2x2 unitary matrix.

This function is 10 times faster than `_random_unitary(2)`.
"""
function _random_unitary_2x2()
    sp_squared = rand() # sin^2(phi)
    sp = sqrt(sp_squared) # sin(phi)
    cp = sqrt(1 - sp_squared) # cos(phi)
    (alpha, beta, gamma) = 2pi .* (rand(), rand(), rand())
    return [cis(beta + alpha)    * cp  cis(gamma + alpha) * sp
            -cis(-gamma + alpha) * sp  cis(-beta + alpha) * cp]
end

function _random_unitary_2x2B()
    sp_squared = rand() # sin^2(phi)
    sp = sqrt(sp_squared) # sin(phi)
    cp = sqrt(1 - sp_squared) # cos(phi)
    (alpha, beta, gamma) = 2 .* (rand(), rand(), rand())
    return [cispi(beta + alpha)    * cp  cispi(gamma + alpha) * sp
            -cispi(-gamma + alpha) * sp  cispi(-beta + alpha) * cp]
end

function _random_unitary_2x2!(a::Matrix{<:Complex})
    size(a) == (2, 2) || throw(DimensionMismatch("Expecting a 2x2 complex matrix"))
    # TODO: check for 1 based, or other layout checks
    sp_squared = rand() # sin^2(phi)
    sp = sqrt(sp_squared) # sin(phi)
    cp = sqrt(1 - sp_squared) # cos(phi)
    (alpha, beta, gamma) = 2pi .* (rand(), rand(), rand())
    # Inbounds seems not to help
    @inbounds a[1,1] = cis(beta + alpha) * cp
    @inbounds a[2,1] = cis(gamma + alpha) * sp
    @inbounds a[1,2] = -cis(-gamma + alpha) * sp
    @inbounds a[2,2] = cis(-beta + alpha) * cp
    return a
end

# There ought to be a way to multiply by the diagonal to save operations.
# But everything I try is more than an order of magnitude slower than
# the naive/simple implementation.
function _random_unitary(n::Integer)
    z = randn(ComplexF64, n, n)
    qr_fac = LinearAlgebra.qr(z)
    ph = LinearAlgebra.diagm([x / abs(x) for x in LinearAlgebra.diag(qr_fac.R)])
    return LinearAlgebra.mul!(z, qr_fac.Q, ph) # reuse storage in z for result
end

function _random_unitary!(A::AbstractMatrix)
    Random.randn!(A)
    qr_fac = LinearAlgebra.qr(A)
    ph = LinearAlgebra.diagm([x / abs(x) for x in LinearAlgebra.diag(qr_fac.R)])
    return LinearAlgebra.mul!(A, qr_fac.Q, ph)
end

# TODO: could make a much faster special case for 2x2 matrices
#
# You can verify with this function that exp(im*theata*U) has det=1
# if U is Hermitian, unitary, with trace 0.
# Obvious for n=2. You can test here for larger n. This must be provable
# This may be good enough if you need some kind of random Hermitian unitary
#  But, this is uniform over eigenvalues. This is probably not
# a good measure. All +1 or all -1 will give +-I, which, in some kind of Haar
# sense, should never happen.
"""
    random_unitary_hermitian(n::Integer; trace=nothing)

Return a random Hermitian, unitary matrix. If `trace` is an
integer, then constrain the trace to this value.

If the trace is not constrained, then you are likely to get the identity operator even
though there is only one such operator and an infinite number of
trace-zero matrices (for even `n`).
"""
function random_unitary_hermitian(n::Integer; trace::Union{Integer,Nothing}=nothing)
    if n < 0
        throw(DomainError(n, "n must be non-negative"))
    end
    m = random_unitary(n)
    # Unitaries have eigvals with modulus 1. Hermitian means real eigvals.
    # So they must be all +-1
    if isnothing(trace)
        d = LinearAlgebra.diagm([rand(Bool) ? -1 : 1 for i in 1:n])
    else
        if ! iseven(trace + n)
            throw(ArgumentError("Requested trace must be of the same parity as matrix dimension n."))
        end
        if trace > n || trace < n
            throw(ArgumentError("Trace must be less than or equal to n and greater than or equal to -n."))
        end
        num_plus_one = div(trace+n, 2)
        d = LinearAlgebra.diagm(Random.shuffle!([i<=num_plus_one ? 1 : -1 for i in 1:n]))
    end
    return m * d * m'
end

# This seems to work, but it's slow
"""
    random_special_unitary(n::Integer)

Return an `n`x`n` random special unitary matrix.

This is computed by scaling the eigenvalues of a Haar-random unitary matrix.
"""
function random_special_unitary(n::Integer)
    u = random_unitary(n)
    eigs = LinearAlgebra.eigen(u)
    (U, evals) = (eigs.vectors, eigs.values)
    scaled_evals = evals ./ prod(evals)^(1/n)
    U' * LinearAlgebra.diagm(scaled_evals) * U
end

"""
    random_normal(n::Integer)

Return a random normal matrix. Eigenvalues have
normally distributed real and imaginary components
"""
function random_normal(n::Integer)
    eigvals = randn(ComplexF64, n)
    U = random_unitary(n)
    U * Diagonal(eigvals) * U'
end
