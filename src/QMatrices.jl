"""
    module QMatrices

Plain-data matrices and vectors for quantum information.

- `I2`: 2q identity
- `H`: Hadamard
- `X, Y, Z`: Pauli matrices
- `SX`: Square root of `X`
- `AB`: A ⊗ B for all pairs of `I,X,Y,Z,H`

Others, not yet organized
```
z0, z1, plus, minus, iplus, iminus,
I2, CH, CCH, S, T,
CX, CCX, CY, CCY, CZ, CCZ,
SWAP, CSWAP, iSWAP, ECR,
Rphi, Rphipi,
RX, RXpi,
RY, RYpi,
RZ, RZpi,
RXX, RXXpi,
RXZ, RXZpi,
RZX, RZXpi,
U, Ualt, U2,
R, Rpi,
magic_basis,
random_unitary,
random_unitary!,
random_special_unitary,
random_unitary_hermitian,
random_normal,
random_orthogonal
```
"""
module QMatrices

export z0, z1, plus, minus, iplus, iminus

export I2, H, CH, CCH, S, T,
    X, CX, CCX, Y, CY, CCY, Z, CZ, CCZ,
    SWAP, CSWAP, iSWAP, ECR,
    Rphi, Rphipi,
    SX,
    RX, RXpi,
    RY, RYpi,
    RZ, RZpi,
    RXX, RXXpi,
    RXZ, RXZpi,
    RZX, RZXpi,
    U, Ualt, U2,
    R, Rpi,
    magic_basis,
    random_unitary,
    random_unitary!,
    random_special_unitary,
    random_unitary_hermitian,
    random_normal,
    random_orthogonal

include("basic_dynamic.jl")
include("linalg.jl")
include("construction.jl")
include("parametric.jl")
include("constructed.jl")
include("random.jl")

# For compiling workflows for statically-compiled-like latency
# using PrecompileTools: @setup_workload, @compile_workload
# include("precompile.jl")

end # module QMatrices
