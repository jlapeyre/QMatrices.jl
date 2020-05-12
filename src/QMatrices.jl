"""
    module QMatrices

Plain-data matrices and vectors for quantum information.
`k0, k1, Id2, X, Y, Z, H, S, T, sqrt_NOT, CX, CCX, CZ, SWAP`
"""
module QMatrices

export k0, k1, Id2, X, Y, Z, H, S, T,
    sqrt_NOT,
    CX, CCX, CY, CZ, SWAP, CSWAP,
    Rphi, Rphipi,
    RX, RXpi,
    RY, RYpi,
    RZ, RZpi

export swap, proj

include("linalg.jl")
include("construction.jl")
include("parametric.jl")
include("basic_dynamic.jl")
include("basic_static.jl")

# using .Dynamic: Rphi, Rphipi, RX, RXpi, RY, RYpi,
#     RZ, RZpi


end # module QMatrices
