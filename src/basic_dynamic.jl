module Dynamic

using QMatrices: checksquare, control, flipbit

export k0, k1, Id2, X, Y, Z, H, S, T,
    CX, CY, CZ, SWAP,
    Rz, Rzpi

# TODO: We should probably define these in a straightforward way, rather than preserve
# generality in the definitions.
"""
    k0

`k0` is `[1, 0]`.
"""
const k0 = [1, 0]

"""
    k1

`k1` is `[0, 1]`.
"""
const k1 = [0, 1]


####
#### Single qubit gates
####

"""
    Id2

The single-qubit identity operator.
"""
const Id2 = [1 0; 0 1]

"""
    X

The Pauli `X` gate.
"""
const X = [0 1; 1 0]

"""
    Y

The Pauli `Y` gate.
"""
const Y = [0 -im; im 0]

"""
    Z

The Pauli `Z` gate.
"""
const Z = [1 0; 0 -1]

"""
    H

The single-qubit Hadamard gate.
"""
const H = (1/sqrt(2)) * [1 1; 1 -1]

"""
    S

The phase gate.
"""
const S = [1 0; 0 im]

"""
    T

The π/8 gate.
"""
const T = [1 0; 0 exp(im * pi / 4)]

# FIXME: These are easy to construct. Maybe we do not need
# to provide them. Decide.
"""
    sqrt_NOT

The square root of the `NOT` (or `X`) gate.
"""
const sqrt_NOT = sqrt(complex(X))

# FIXME: Highlighting of $cgate is broken
# Controlled X, Y, Z with one and two control qubits
for gate in (:X, :Y, :Z, :H)
    for n in 1:2
        cgate = Symbol(fill(:C, n)..., gate)
        gatestr = string(gate)
        cgatestr = string(cgate)
        expl = n == 2 ? "gate with two control qubits." : "gate."
        @eval begin
"""
   $($cgatestr)

The controlled-`$($gatestr)` $($expl)
"""
        const $cgate = control($gate, $n)
        end
    end
end

"""
    SWAP


The 2-qubit `SWAP` gate.
"""
const SWAP = [1 0 0 0
              0 0 1 0
              0 1 0 0
              0 0 0 1]

"""
    CSWAP

The controlled `SWAP` gate.
"""
const CSWAP = control(SWAP)

end # module Dynamic
