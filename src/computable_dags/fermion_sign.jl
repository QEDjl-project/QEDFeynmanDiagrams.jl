"""
    relative_sign_pair(p1::VirtualParticle, p2::VirtualParticle)

Returns true if the Pair task combining p1 and p2 into res needs a relative sign in the DAG.
"""
function relative_sign_pair(p1::VirtualParticle, p2::VirtualParticle)
    if (p1.species == Photon || p2.species == Photon)
        #@info "P [FALS] relative sign of $p1 + $p2 -> $res"
        return false
    end

    local vec::Vector{OPEN_FERMION_CYCLE_T}

    (l, el) = if p1.species == Electron
        _canonical_index(p1)
    elseif p2.species == Electron
        _canonical_index(p2)
    end
    (r, po) = if p1.species == Positron
        _canonical_index(p1)
    elseif p2.species == Positron
        _canonical_index(p2)
    end

    vec = [p1.open_cycles..., p2.open_cycles..., (el, po)]

    n = _count_closed_cycles(vec)

    return n % 2 == 1
end

"""
    relative_sign_triple(p2::VirtualParticle, p3::VirtualParticle, p3::VirtualParticle)

Returns true if the Triple task combining electron p1, and positron p2 into a full diagram family needs a relative sign in the DAG.
"""
function relative_sign_triple(p1::VirtualParticle, p2::VirtualParticle, p3::VirtualParticle)
    # p1 + p2 is the photon
    n1 = relative_sign_pair(p1, p2)

    new_cycle = (_canonical_index(p1)[2], _canonical_index(p2)[2])

    res = VirtualParticle(
        p1.proc,
        Photon(),
        (_contributions(p1) + _contributions(p2))...,
        sort(
            reduce_cycles(
                OPEN_FERMION_CYCLE_T[p1.open_cycles..., p2.open_cycles..., new_cycle]
            ),
        ),
    )

    # cycles closed?
    closed_cycle_count = _count_closed_cycles(
        OPEN_FERMION_CYCLE_T[p3.open_cycles..., res.open_cycles...]
    )

    n2 = closed_cycle_count % 2 == 1

    return n1 != n2 # xor
end
