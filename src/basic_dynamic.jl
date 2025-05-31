const z0 = [1, 0]
const z1 = [0, 1]
const plus = (z0 + z1) / sqrt(2)
const minus = (z0 - z1) / sqrt(2)
const iplus = (z0 + im * z1) / sqrt(2)
const iminus = (z0 - im * z1) / sqrt(2)

"""
    z0, z1, plus, minux, iplus, iminus

The states |0>, |1>, |+⟩, |-⟩, |+i⟩, |-i⟩. These are eigenstates of the Pauli Z, X, and Y operators.

- `z0` is the `+1` eigenvector, of the Pauli [`Z`](@ref) operator.
- `z1` is the `-1` eigenvector, of the Pauli [`Z`](@ref) operator.
- `minus` is the `-1` eigenvector, of the Pauli [`X`](@ref) operator.
- `plus` is the `+1` eigenvector, of the Pauli [`X`](@ref) operator.
- `iminus` is the `-1` eigenvector, of the Pauli [`Y`](@ref) operator.
- `iplus` is the `+1` eigenvector, of the Pauli [`Y`](@ref) operator.

# Examples

```jldoc-test
julia> X * (0.1 * plus + 0.2 * minus) == (0.1 * plus - 0.2 * minus)
true

julia> Y * (0.1 * iplus + 0.2 * iminus) == (0.1 * iplus - 0.2 * iminus)
true

julia> Z * (0.1 * z0 + 0.2 * z1) == (0.1 * z0 - 0.2 * z1)
true
```
"""
z0, z1, plus, minus, iplus, iminus

####
#### Single qubit gates
####

"""
    I2

The single-qubit identity operator.
"""
const I2 = [1 0; 0 1]

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
const T = [1 0; 0 1/sqrt(2) * (1 + im)]
# Following slightly less accurate.
# Or better said; above T^2 - S == 0. But not for following
# const T = [1 0; 0 (cospi(1/4) + im * sinpi(1/4))]

# In recent versions of Julia, computing this at compile time takes several
# seconds. So we insert the result.
"""
    sqrt_NOT

The square root of the `NOT` (or `X`) gate.
"""
const sqrt_NOT = ComplexF64[1.0 + 1.0im 1.0 - 1.0im; 1.0 - 1.0im 1.0 + 1.0im]
# = sqrt(complex(X))

const SX = 1//2 * [1+im 1-im; 1-im 1+im]


"""
    SWAP


The 2-qubit `SWAP` gate.
"""
const SWAP = [1 0 0 0
              0 0 1 0
              0 1 0 0
              0 0 0 1]

"""
    iSWAP


The 2-qubit `iSWAP` gate.
"""
const iSWAP = [1 0  0  0
               0 0  im 0
               0 im 0  0
               0 0  0  1]

"""
    ECR

Echoed cross-resonance gate.
"""
const ECR = (1/sqrt(2)) * (kron(I2, X) - kron(X, Y))

"""
    Magic basis

The matrix whose columns are the magic basis.
"""
const magic_basis = (1 / sqrt(2)) *
    [1  0    0   im
     0  im   1    0
     0  im  -1    0
     1  0    0  -im]
