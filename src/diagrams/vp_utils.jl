"""
    disjunct(a::VirtualParticle, b::VirtualParticle)

Return true if the momenta contributions of `a` and `b` are disjunct.
"""
function disjunct(a::VirtualParticle, b::VirtualParticle)
    for (a_contrib, b_contrib) in
        Iterators.zip(Iterators.flatten.(_contributions.((a, b)))...)
        if b_contrib && a_contrib
            return false
        end
    end

    return true
end

"""
    contains(a::VirtualParticle, b::VirtualParticle)

Returns true if the set of particles contributing to `a` are contains the set of particles contributing to `b`.
"""
function contains(a::VirtualParticle, b::VirtualParticle)
    for (a_contrib, b_contrib) in
        Iterators.zip(Iterators.flatten.(_contributions.((a, b)))...)
        if b_contrib && !a_contrib
            return false
        end
    end

    return true
end

"""
    are_total(a::VirtualParticle, b::VirtualParticle, c::VirtualParticle)

Return true if a, b and c combined contain all external particles exactly once.
"""
function are_total(
    a::VirtualParticle{PROC}, b::VirtualParticle{PROC}, c::VirtualParticle{PROC}
) where {PROC<:AbstractProcessDefinition}
    for (a_contrib, b_contrib, c_contrib) in
        Iterators.zip(Iterators.flatten.(_contributions.((a, b, c)))...)
        if a_contrib + b_contrib + c_contrib != 1
            return false
        end
    end

    (l, el) = if a.species == Electron
        _canonical_index(a)
    elseif b.species == Electron
        _canonical_index(b)
    else
        _canonical_index(c)
    end
    (r, po) = if a.species == Positron
        _canonical_index(a)
    elseif b.species == Positron
        _canonical_index(b)
    else
        _canonical_index(c)
    end

    reduced_cycles = reduce_cycles(
        OPEN_FERMION_CYCLE_T[a.open_cycles..., b.open_cycles..., c.open_cycles..., (el, po)]
    )
    # if the combination is total, there cannot be any leftover open cycles
    if reduced_cycles != OPEN_FERMION_CYCLE_T[]
        @warn "rejected!"
        return false
    end

    return true
end

@inline _make_up_helper(a::Tuple{}, b::Tuple{}, c::Tuple{}) = true
@inline function _make_up_helper(
    a::Tuple{Bool,Vararg}, b::Tuple{Bool,Vararg}, c::Tuple{Bool,Vararg}
)
    return if a[begin] + b[begin] == c[begin]
        _make_up_helper(a[2:end], b[2:end], c[2:end])
    else
        false
    end
end

"""
    make_up(a::VirtualParticle, b::VirtualParticle, c::VirtualParticle)
    
For virtual particles `a`, `b`, and `c`, return true if `a` and `b`'s joint momentum contributions add up to `c`'s momentum contributions.
"""
function make_up(
    a::VirtualParticle{PROC,I,O}, b::VirtualParticle{PROC,I,O}, c::VirtualParticle{PROC,I,O}
) where {PROC,I,O}
    if a.species == b.species == Photon
        return false
    end

    if !_make_up_helper(_in_contributions(a), _in_contributions(b), _in_contributions(c))
        return false
    end
    if !_make_up_helper(_out_contributions(a), _out_contributions(b), _out_contributions(c))
        return false
    end

    cycles = if a.species != Photon && b.species != Photon
        (s1, n1) = _canonical_index(a)
        (s2, n2) = _canonical_index(b)

        new_cycle = if (s1 == :left && s2 == :right)
            (n1, n2)
        elseif (s1 == :right && s2 == :left)
            (n2, n1)
        else
            @assert false
        end

        reduce_cycles(OPEN_FERMION_CYCLE_T[a.open_cycles..., b.open_cycles..., new_cycle])
    else
        reduce_cycles(OPEN_FERMION_CYCLE_T[a.open_cycles..., b.open_cycles...])
    end

    if sort(cycles) != c.open_cycles
        return false
    end

    return true
end
