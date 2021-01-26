export E
export b00, b01, b10, b11, plus, minus, iplus, iminus

"""
    E(n::Integer, m::Integer=n, a::Integer, b::Integer)

Computational-basis vectors of `n`x`m` matrices.
"""
function E(n::Integer, m::Integer, a::Integer, b::Integer)
    (1 <= a <= n && 1 <= b <= m) || error("E: Invalid indices, violate 1 <= $a <= $n && 1 <= $b <= $m")
    Em = zeros(Int8,  n, m)
    Em[a, b] = 1
    return Em
end
E(n::Integer, a::Integer, b::Integer) = E(n, n, a, b)

####
#### Bell States
####

"""
    b00, b01, b10, b11

Bell states.
"""
const b00 = (ket(0, 0) + ket(1, 1)) / sqrt(2)
const b10 = (ket(0, 0) - ket(1, 1)) / sqrt(2)
const b01 = (ket(0, 1) + ket(1, 0)) / sqrt(2)
const b11 = (ket(0, 1) - ket(1, 0)) / sqrt(2)
@doc (@doc b00) b01
@doc (@doc b00) b10
@doc (@doc b00) b11

####
#### Plus and minus states
####

"""
    plus, minux, iplus, iminus

The states |+⟩, |-⟩, |+i⟩, |-i⟩. These are eigenstates of
the Pauli X and Y operators.
"""
const plus = (ket(0) + ket(1)) / sqrt(2)
const minus = (ket(0) - ket(1)) / sqrt(2)
const iplus = (ket(0) + im * ket(1)) / sqrt(2)
const iminus = (ket(0) - im * ket(1)) / sqrt(2)
@doc (@doc plus) minus
@doc (@doc plus) iplus
@doc (@doc plus) iminus
