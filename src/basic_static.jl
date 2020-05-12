using StaticArrays: SArray

# We define the corresponding heap-allocated quantities first, then convert them
# to these static quantities.
# This is, in part, because defining the static quantities first requires a lot of,
# and varied, boilerplate.

"""
    sarray(A)

Convert `A` to a `StaticArrays.SArray`.
"""
sarray(A) = SArray{Tuple{size(A)...}}(A)

# Create static versions of heap-allocated objects
for A in (:k0, :k1, :X, :Y, :Z, :Id2, :H,
          :S, :T, :sqrt_NOT, :CX, :CCX, :CY, :CZ, :SWAP, :CSWAP)
    @eval const $A = sarray(Dynamic.$A)
    @eval @doc (@doc Dynamic.$A) $A
end
