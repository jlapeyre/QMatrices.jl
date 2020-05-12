# These are copied from LinarAlgebra because they are not exported,
# and are perhaps not part of the API. This should be checked.
# Julia's API definition is semi-formal and not easily verified.

"""
    QMatrices.checksquare(A)

Check that a matrix is square, then return its common dimension.
For multiple arguments, return a vector.

# Examples
```jldoctest
julia> A = fill(1, (4,4)); B = fill(1, (5,5));

julia> LinearAlgebra.checksquare(A, B)
2-element Array{Int64,1}:
 4
 5
```
"""
function checksquare(A)
    m,n = size(A)
    m == n || throw(DimensionMismatch("matrix is not square: dimensions are $(size(A))"))
    m
end

function checksquare(A...)
    sizes = Int[]
    for a in A
        size(a,1)==size(a,2) || throw(DimensionMismatch("matrix is not square: dimensions are $(size(a))"))
        push!(sizes, size(a,1))
    end
    return sizes
end
