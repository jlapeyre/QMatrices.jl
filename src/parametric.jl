####
#### Parametric gates
####

"""
    Rphi(ϕ)

The phase shift gate ``R_ϕ``.
"""
function Rphi(ϕ)
    return [1 0; 0 exp(im * ϕ)]
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
    return [1 0; 0 _exp_ipi(z)]
end

function RX(θ)
    c = cos(θ/2)
    s = - im *sin(θ/2)
    return [c s; s c]
end

function RXpi(z)
    c = cospi(z/2)
    s = -im * sinpi(z/2)
    return [c s; s c]
end

function RY(θ)
    c = cos(θ/2)
    s = sin(θ/2)
    return [c -s; s c]
end

function RYpi(z)
    c = cospi(z/2)
    s = sinpi(z/2)
    return [c -s; s c]
end

function RZ(θ)
    return [exp(-im*θ/2) 0; 0 exp(im*θ/2)]
end

function RZpi(z)
    return [_exp_ipi(-z/2) 0; 0 _exp_ipi(z/2)]
end

"""
    U(θ, ϕ, λ)

Matrix from SU(2).
"""
function U(θ, ϕ, λ)
    c = cos(θ/2)
    s = sin(θ/2)
    fpl =  (ϕ + λ) / 2
    fml =  (ϕ - λ) / 2
    cfpl = cos(fpl)
    sfpl = sin(fpl)
    cfml = cos(fml)
    sfml = sin(fml)
    f00 = complex(cfpl, -sfpl) * c
    f01 = complex(-cfml, sfml) * s
    f10 = complex(cfml, sfml) * s
    f11 = complex(cfpl, sfpl) * c
    return [f00 f01; f10 f11]
end

"""
    Upi(θ, ϕ, λ)

Matrix from SU(2), with `θ`, `ϕ`, and `λ` given as multiples of `π`.
This is more accurate than `U`.
"""
function Upi(θ, ϕ, λ)
    c = cospi(θ/2)
    s = sinpi(θ/2)
    fpl =  (ϕ + λ) / 2
    fml =  (ϕ - λ) / 2
    cfpl = cospi(fpl)
    sfpl = sinpi(fpl)
    cfml = cospi(fml)
    sfml = sinpi(fml)
    f00 = complex(cfpl, -sfpl) * c
    f01 = complex(-cfml, sfml) * s
    f10 = complex(cfml, sfml) * s
    f11 = complex(cfpl, sfpl) * c
    return [f00 f01; f10 f11]
end
