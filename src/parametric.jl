####
#### Parametric gates
####

# TODO!!
# All the methods that construct a matrix with elements in place are insanely slow.
# Instead allocate u matrix with undef elements and set them one by one.
# This is also quite slow (not as slow) unless you explicitly set them one at a time with a literal
# index. This is very unfortunate, but I don't know a way around it at the moment.

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

function _R(θ, ϕ, sincosf)
    (s, c) = sincosf(θ/2)
    (sp, cp) = sincosf(ϕ)

    m = Matrix{Complex{typeof(θ)}}(undef, 2, 2)
    @inbounds m[1] = c
    @inbounds m[2] = (-im*cp + sp) * s
    @inbounds m[3] = (-im*cp - sp) * s
    @inbounds m[4] = c
    return m

    # Slower :(
    # for (i, x) in enumerate([c, (-im*cp + sp) * s, (-im*cp - sp) * s, c])
    #     @inbounds m[i] = x
    # end
    # m
    # More slower. Following is several times slower method in use above
    # return [c  (-im*cp - sp) * s;
    #         (-im*cp + sp) * s c]
end

@doc raw"""
    R(θ, ϕ)

The R gate.
```math
R(\theta, \phi) = e^{-i \frac{\theta}{2} (\cos{\phi} x + \sin{\phi} y)}
```
"""
function R(θ, ϕ)
    _R(θ, ϕ, sincos)
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
    _R(θ, ϕ, sincospi)
end

"""
    RXX(θ)

The R_XX gate, `exp(-im θ / 2 X⊗X)`.
"""
function RXX(θ)
    _RXX(θ, sincos)
end

"""
    RXXpi(θ)

The RXX gate scaled by π.
`exp(-im θ /(2π) X⊗X)`.
"""
function RXXpi(θ)
    _RXX(θ, sincospi)
end

function _RXX(θ, sincosf)
    (s_a, c) = sincosf(θ / 2)
    s = -im * s_a
    z = zero(θ)
    [c z z s
     z c s z
     z s c z
     s z z c]
end

"""
    RXXYY(θ)

The R_XXYY gate, `exp(-im θ / 2 (X⊗X + Y⊗Y))`.
"""
function RXXYY(θ)
    _RXXYY(θ, cos, sin)
end

"""
    RXXYYpi(θ)

The R_XXYY gate with `θ` scaled by `π`. `exp(-im * pi * θ / 2 * (X⊗X + Y⊗Y))`
"""
function RXXYYpi(θ)
    _RXXYY(θ, cospi, sinpi)
end


function _RXXYY end

let
    IIh = II / 2
    ZZh = ZZ / 2
    IIZZp = IIh + ZZh
    IIZZm = IIh - ZZh
    XXYY = (XX + YY)/2
    global _RXXYY
    function _RXXYY(θ, cfunc, sfunc)
        return IIZZp + cfunc(θ) * IIZZm -im * sfunc(θ) * XXYY
    end
end

"""
    RXZ(θ)

`exp(-im * θ / 2 X ⊗ Z)`
"""
function RXZ(θ)
    _RXZ(θ, cos, sin)
end

"""
    RXZ(θ)

`exp(-im * θ / (2 π) X ⊗ Z)`
"""
function RXZpi(θ)
    _RXZ(θ, cospi, sinpi)
end

function _RXZ(θ, cosf, sinf)
    t =  θ / 2
    c = cosf(t)
    s = im * sinf(t)
    [c  0 -s  0
     0  c  0  s
    -s  0  c  0
     0  s  0  c]
end

"""
    RZZ(θ)

`exp(-im * θ/2 Z ⊗ Z)`
"""
function RZZ(θ)
    _RZZ(θ, sincos)
end

"""
    RZZpi(θ)

`exp(-im * pi * θ/2 Z ⊗ Z)`
"""
function RZZpi(θ)
    _RZZ(θ, sincospi)
end

# This is much faster than writing out the elements explicitly
# I don't understand why.
# There are far fewer alloctations this way, compared to writing
# them explicitly. That should actually only allocate once.
function _RZZ(θ, sincosf)
    t =  θ / 2
    (s, c) = sincosf(t)
    c * II - im * s * ZZ
end

# In Qiskit, there is only RXZ, but they write
# two versions with swapped qubits.
# TODO: Make sure endianness story is straight here.
"""
    RZX(θ)

`exp(-im * θ / 2 Z ⊗ X)`
"""
function RZX(θ)
    _RZX(θ, cos, sin)
end

"""
    RZX(θ)

`exp(-im * θ / (2 π) X ⊗ Z)`
"""
function RZXpi(θ)
    _RZX(θ, cospi, sinpi)
end

function _RZX(θ, cosf, sinf)
    t =  θ / 2
    c = cosf(t)
    s = im * sinf(t)
    [c -s  0  0
    -s  c  0  0
     0  0  c  s
     0  0  s  c]
end
