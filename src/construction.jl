import LinearAlgebra
import ILog2

export flipbit, control
export swap, projector, ket, ketbra, zeroket, oneket
export ⊗

const ⊗ = kron

# FIXME: I think I copied this from somewhere. Probably get rid of it.
"""
    flipbit(x)

Return `1` if `x` is `0` and vice versa.
"""
flipbit(x) = x == 1 ? zero(x) : one(x)

"""
    swap(G::AbstractMatrix)

Swap the first and second qubits in the two-qubit operator `G`.
"""
swap(G) = SWAP * G * SWAP

"""
    ket(b::Integer)
    ket(bs::Integer...)
    ket(s::AbstractString)

Return a vector representation the single-qubit computational basis state indexed by `b`, either `0` or `1`.
This returns either `z0` or `z1`.

Return the `length(bs)`-qubit computational basis state as a one-hot vector. Each argument must be either `0` or `1`.

Construct a computational-basis ket from bitstring `s`.
"""
ket(b::Integer) = b == 0 ? z0 : b == 1 ? z1 : throw(DomainError(b))

function _ket(position::Integer, n_qubits::Integer)
    v = zeros(Int, 2^n_qubits)
    v[position] = 1
    return v
end

ket(b0::Integer, bs::Integer...) = _ket_from_int_args(b0, bs...)

function _ket_from_int_args(bs::Integer...)
    pos = 1
    n = length(bs)
    for (i, val) in enumerate(bs)
        pos += val << (n-i)
    end
    return _ket(pos, length(bs))
end

# Example: ket("01101")
function ket(s::AbstractString)
    pos = parse(Int, s, base=2) + 1
    n_qubits = length(s)
    return _ket(pos, n_qubits)
end

"""
    zeroket(n_qubits)

Return the `n_qubits`-qubit zero ket.
"""
zeroket(n_qubits::Integer) = _ket(1, n_qubits)

"""
    oneket(n_qubits)

Return the `n_qubits`-qubit one ket.
"""
oneket(n_qubits::Integer) = _ket(2^n_qubits, n_qubits)

"""
    ketbra(ket::AbstractVector, bra::LinearAlgebra.Adjoint)

Return `ket * bra`. That is the matrix element `|ket><bra|`.
"""
ketbra(ket::AbstractVector, bra::LinearAlgebra.Adjoint) = ket * bra

"""
    ketbra(ket1::AbstractVector, ket2::AbstractVector)

Return `ket1 * ket2'`
"""
ketbra(ket1::AbstractVector, ket2::AbstractVector) = ketbra(ket1, ket2')

# TODO: This could be hardcoded rather than calculated
"""
    ketbra(b0::Integer, b1::Integer)

Return `ketbra(ket(b0), ket(b1))`
"""
ketbra(b0::Integer, b1::Integer) = ketbra(ket(b0), ket(b1))

# Use a materialized covector. This is 8 or so times slower than using the adjoint wrapper (On Julia 1.2)
ketbra(ket1::AbstractVector, ket2::AbstractMatrix) = ket1 * ket2 # Assume dimensions are correct. i.e. ket2 is an adjoint

"""
    projector(ket::AbstractVector)

Return `ketbra(ket, ket)`.
"""
projector(ket::AbstractVector) = ketbra(ket, ket)

"""
    projector(ket::Integer)

Return `projector(ket(b))`.
"""
projector(b::Integer) = projector(ket(b))

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

"""
    E(n::Integer, m::Integer=n, a::Integer, b::Integer)

Computational-basis vectors of `n`x`m` matrices.
"""
function E(n::Integer, m::Integer, a::Integer, b::Integer)
    (1 <= a <= n && 1 <= b <= m) || error("E: Invalid indices, violate 1 <= $a <= $n && 1 <= $b <= $m")
    Em = zeros(Int8,  m, m)
    Em[a, b] = 1
    return Em
end
E(n::Integer, a::Integer, b::Integer) = E(n, n, a, b)


####
#### Bell States
####

const b00 = (ket(0, 0) + ket(1, 1)) / sqrt(2)
const b10 = (ket(0, 0) - ket(1, 1)) / sqrt(2)
const b01 = (ket(0, 1) + ket(1, 0)) / sqrt(2)
const b11 = (ket(0, 1) - ket(1, 0)) / sqrt(2)
