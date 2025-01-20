export E
export b00, b01, b10, b11, plus, minus, iplus, iminus
export bell_proj, bell, bell_diag, comp_to_bell, bell_to_comp

"""
    E(n::Integer, m::Integer=n, a::Integer, b::Integer)

Computational-basis vectors of `n`x`m` matrices.
"""
function E(n::Integer, m::Integer, a::Integer, b::Integer)
    (1 <= a <= n && 1 <= b <= m) || error("E: Invalid indices, violate 1 <= $a <= $n && 1 <= $b <= $m")
    Em = zeros(Int8,  n, m)
    Em[a, b] = 1
    return Em
end
E(n::Integer, a::Integer, b::Integer) = E(n, n, a, b)

####
#### Bell States
####

"""
    b00, b01, b10, b11

Bell states.
"""
const b00 = (ket(0, 0) + ket(1, 1)) / sqrt(2)
const b10 = (ket(0, 0) - ket(1, 1)) / sqrt(2)
const b01 = (ket(0, 1) + ket(1, 0)) / sqrt(2)
const b11 = (ket(0, 1) - ket(1, 0)) / sqrt(2)
@doc (@doc b00) b01
@doc (@doc b00) b10
@doc (@doc b00) b11

# indexable Bell states
const _bells = (b00, b01, b10, b11)

"""
    bell_proj

Tuple of four projectors onto Bell states.
"""
const bell_proj = proj.(_bells)

"""
    bell(a, b)

Return the Bell state ``b_{a,b}``. `a` and `b` take values `0` and `1`.
"""
function bell(a, b)
    return _bells[2*a + b + 1]
end


"""
    bell_diag(c0, c1, c2, c3)

Return a Bell-diagonal state in the computational basis with
the standard encoding.
"""
bell_diag(c0, c1, c2, c3) = sum((c0, c1, c2, c3) .* bell_proj)


"""
    comp_to_bell

Unitary transform from computational to Bell basis.
"""
const comp_to_bell = sum(bell(a,b) * ket(a,b)' for a in (0,1) for b in (0,1))

"""
    bell_to_comp

Unitary transform from Bell to computational basis.
"""
const bell_to_comp = collect(comp_to_bell')


# FIXME: Highlighting of $cgate is broken
# Controlled X, Y, Z, H with one and two control qubits
for gate in (:X, :Y, :Z, :H)
    for n in 1:2
        cgate = Symbol(fill(:C, n)..., gate)
        gatestr = string(gate)
        cgatestr = string(cgate)
        expl = n == 2 ? "gate with two control qubits." : "gate."
        @eval begin
"""
   $($cgatestr)

The controlled-`$($gatestr)` $($expl)
"""
        const $cgate = control($gate, $n)
        end
    end
end


"""
    CSWAP

The controlled `SWAP` gate.
"""
const CSWAP = control(SWAP)

####
#### Plus and minus states
####

