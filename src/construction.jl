import ILog2

export flipbit, control
export ⊗

const ⊗ = kron

"""
    flipbit(x)

Return `1` if `x` is `0` and vice versa.
"""
flipbit(x) = x == 1 ? zero(x) : one(x)

swap(G) = SWAP * G * SWAP

proj(v::AbstractVector) = v * v'

"""
    control(A::AbstractMatrix, n_control_bits=1)

Return a controlled-`A` gate with `n_control_bits` control bits.
"""
function control(A::AbstractMatrix, n_control_bits=1)
    T = eltype(A)
    dA = checksquare(A)
    nqubits = ILog2.ilog2(dA) + n_control_bits
    d = 2^nqubits
    controlled_A = zeros(T, (d, d))
    @inbounds for i in 1:(d-dA)
        controlled_A[i, i] = one(T)
    end
    @inbounds controlled_A[d-dA+1:d, d-dA+1:d] .= A
    return controlled_A
end
