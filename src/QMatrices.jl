"""
    module QMatrices

Plain-data matrices and vectors for quantum information.
`z0, z1, Id2, X, Y, Z, H, S, T, sqrt_NOT, CX, CCX, CZ, SWAP`
"""
module QMatrices

export z0, z1, Id2, H, CH, CCH, S, T,
    sqrt_NOT,
    X, CX, CCX, Y, CY, CCY, Z, CZ, CCZ,
    SWAP, CSWAP,
    Rphi, Rphipi,
    RX, RXpi,
    RY, RYpi,
    RZ, RZpi,
    U, Upi,
    E,
    b00, b01, b10, b11

include("linalg.jl")
include("construction.jl")
include("parametric.jl")
include("basic_dynamic.jl")

end # module QMatrices
