####
#### Parametric gates
####

"""
    Rphi(ϕ)

The phase shift gate ``R_ϕ``.
This is equivalent to qiskit's `u1` gate.
"""
function Rphi(ϕ)
    return [1 0; 0 cis(ϕ)]
end

# Could use complex(reverse(sincospi(z))...)
"""
    _exp_ipi(z)

Return `exp(im * pi * z)`. This uses the accurate `cospi` and `sinpi` functions.
"""
_exp_ipi(z) = cispi(z)

"""
    Rphipi(z)

Return the ``R_ϕ`` gate for ``ϕ = zπ``. This is more accurate
than `Rphi(pi * z)`.
"""
function Rphipi(z)
    return [1 0; 0 _exp_ipi(z)]
end

"""
    U2(ϕ, λ)

u2 gate. Need to put math definitions in here.
"""
function U2(ϕ, λ)
    return [1 -cis(λ); cis(ϕ) cis(ϕ + λ)] / sqrt(2)
end

# This does not appar any faster than computing each separately
# function RXalt(θ)
#     (s, c) = sincos(θ/2)
#     si = - im * s
#     return [c si; si c]
# end

function RX(θ)
    c = cos(θ/2)
    s = - im * sin(θ/2)
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
    return [cis(-θ/2) 0; 0 cis(θ/2)]
end

function RZpi(z)
    return [_exp_ipi(-z/2) 0; 0 _exp_ipi(z/2)]
end

"""
    Ualt(θ, ϕ, λ)

Matrix from SU(2).
Alternative parameterization.
"""
function Ualt(θ, ϕ, λ)
    c = cos(θ/2)
    s = sin(θ/2)
    fpl =  (ϕ + λ) / 2
    cfpl = cos(fpl)
    sfpl = sin(fpl)
    fml =  (ϕ - λ) / 2
    cfml = cos(fml)
    sfml = sin(fml)
    f00 = complex(cfpl, -sfpl) * c
    f01 = complex(-cfml, sfml) * s
    f10 = complex(cfml, sfml) * s
    f11 = complex(cfpl, sfpl) * c
    return [f00 f01; f10 f11]
end

"""
    U(θ, ϕ, λ)

Matrix from SU(2). This is
the same as qiskit's U or U3.
"""
function U(θ, ϕ, λ)
    c = cos(θ/2)
    s = sin(θ/2)
    fpl =  (ϕ + λ)
    cfpl = cos(fpl)
    sfpl = sin(fpl)
    f00 = c
    f01 = -s * complex(cos(λ), sin(λ))
    f10 = s * complex(cos(ϕ), sin(ϕ))
    f11 = c * complex(cfpl, sfpl)
    return [f00 f01; f10 f11]
end

"""
    Ualtpi(θ, ϕ, λ)

Matrix from SU(2), with `θ`, `ϕ`, and `λ` given as multiples of `π`.
This is more accurate than `Ualt`.
"""
function Ualtpi(θ, ϕ, λ)
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

@doc raw"""
    R(θ, ϕ)

The R gate.
```math
R(\theta, \phi) = e^{-i \frac{\theta}{2} (\cos{\phi} x + \sin{\phi} y)}
```
"""
function R(θ, ϕ)
    c = cos(θ/2)
    s = sin(θ/2)
    return [c -im*exp(-im*ϕ)*s;
            -im*exp(im*ϕ)*s c]
end

@doc raw"""
    Rpi(θ, ϕ)

The R gate with θ and ϕ reduced by π. This is may be more accurate,
for example with integral and half-integral multiples of π.
```math
Rpi(\theta, \phi) = e^{-i \frac{\pi\theta}{2} (\cos{\pi\phi} x + \sin{\pi\phi} y)}
```
"""
function Rpi(θ, ϕ)
    c = cospi(θ/2)
    s = sinpi(θ/2)
    return [c -im*cispi(-ϕ)*s;
            -im*cispi(ϕ)*s c]
end
