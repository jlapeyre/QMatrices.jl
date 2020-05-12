####
#### Parametric gates
####

using StaticArrays: @SMatrix

"""
    Rphi(ϕ)

The phase shift gate ``R_ϕ``.
"""
function Rphi(ϕ)
    return @SMatrix [1 0; 0 exp(im * ϕ)]
end

"""
    _exp_ipi(z)

Return `exp(im * pi * z)`. This uses the accurate `cospi` and `sinpi` functions.
"""
_exp_ipi(z) = complex(cospi(z), sinpi(z))

"""
    Rphipi(z)

Return the ``R_ϕ`` gate for ``ϕ = zπ``. This is more accurate
than `Rphi(pi * z)`.
"""
function Rphipi(z)
    return @SMatrix [1 0; 0 _exp_ipi(z)]
end

function RX(θ)
    c = cos(θ/2)
    s = - im *sin(θ/2)
    return @SMatrix [c s; s c]
end

function RXpi(z)
    c = cospi(z/2)
    s = -im * sinpi(z/2)
    return @SMatrix [c s; s c]
end

function RY(θ)
    c = cos(θ/2)
    s = sin(θ/2)
    return @SMatrix [c -s; s c]
end

function RYpi(z)
    c = cospi(z/2)
    s = sinpi(z/2)
    return @SMatrix [c -s; s c]
end

function RZ(θ)
    return @SMatrix [exp(-im*θ/2) 0; 0 exp(im*θ/2)]
end

function RZpi(z)
    return @SMatrix [_exp_ipi(-z/2) 0; 0 _exp_ipi(z/2)]
end
