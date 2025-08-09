const FLOAT_TYPE = Float64

const GATE_NAME_MAP = Dict{Symbol, Matrix{Complex{Float64}}}()

function mybigfloat(z::Complex)
    complex(Float64(z.re), Float64(z.im))
end

function mybigfloat(x)
    Float64(x)
end

for g in (:X, :Y, :Z, :H, :S, :T, :SX)
    m1 = complex(@eval $g)
    gm = mybigfloat.(m1)
    GATE_NAME_MAP[g] = gm
end

GATE_NAME_MAP[:W] = mybigfloat.(cispi(1/4) * [1 0; 0 1])

function gate_from_name(s)
    return GATE_NAME_MAP[Symbol(s)]
end


"""
    compose(gates::AbstractString; map=gate_from_name)

Compose the gates in `gates` from right to left.

`map` is a function mapping a character (or Symbol) to
a dense matrix representation of the gate.

### Examples:
```julia
compose("XYZ")
```
"""
function compose(gates::AbstractString; map=gate_from_name)
    c = complex(float([1 0; 0 1]))
    for g in Iterators.reverse(gates)
        c = gate_from_name(g) * c
    end
    return c
end
