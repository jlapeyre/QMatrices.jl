"""
    module QMatrices

Plain-data matrices and vectors for quantum information.
`z0, z1, I2, X, Y, Z, H, S, T, sqrt_NOT, CX, CCX, CZ, SWAP`
"""
module QMatrices

# For compiling workflows for statically-compiled-like latency
using PrecompileTools: @setup_workload, @compile_workload

export z0, z1, plus, minus, iplus, iminus

export I2, H, CH, CCH, S, T,
    sqrt_NOT,
    X, CX, CCX, Y, CY, CCY, Z, CZ, CCZ,
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

export SX

include("basic_dynamic.jl")
include("linalg.jl")
include("construction.jl")
include("parametric.jl")
include("constructed.jl")
include("random.jl")

include("precompile.jl")

end # module QMatrices
