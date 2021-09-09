"""
    module QMatrices

Plain-data matrices and vectors for quantum information.
`z0, z1, I2, X, Y, Z, H, S, T, sqrt_NOT, CX, CCX, CZ, SWAP`
"""
module QMatrices

export z0, z1, I2, H, CH, CCH, S, T,
    sqrt_NOT,
    X, CX, CCX, Y, CY, CCY, Z, CZ, CCZ,
    SWAP, CSWAP, iSWAP,
    Rphi, Rphipi,
    RX, RXpi,
    RY, RYpi,
    RZ, RZpi,
    U, Upi, Ualt, U2,
    R, Rpi

include("linalg.jl")
include("construction.jl")
include("parametric.jl")
include("basic_dynamic.jl")
include("constructed.jl")

end # module QMatrices
