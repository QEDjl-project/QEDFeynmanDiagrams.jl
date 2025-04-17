"""
    _get_canonical_index(proc::AbstractProcessDefinition, dir::ParticleDirection, index::Int)

Returns a tuple of a symbol which is either `:left` or `:right` in the case of fermions, and `:boson`, in the case of a boson, and an `Int` giving the index of the particle.
"""
function _get_canonical_index(
        proc::AbstractProcessDefinition, dir::ParticleDirection, index::Int
    )
    species = particles(proc, dir)[index]
    species_index = _species_index(proc, dir, species, index)

    inc_parts = number_particles(proc, Incoming(), Electron())
    inc_anti_parts = number_particles(proc, Incoming(), Positron())

    if is_boson(species)
        return (:boson, species_index)
    elseif is_particle(species) && is_incoming(dir)
        return (:left, species_index)
    elseif is_anti_particle(species) && is_outgoing(dir)
        return (:left, species_index + inc_parts)
    elseif is_anti_particle(species) && is_incoming(dir)
        return (:right, species_index)
    elseif is_particle(species) && is_outgoing(dir)
        return (:right, species_index + inc_anti_parts)
    end

    throw("unknown species/dir combination encountered: $(species)/$(dir)")
end

function _species_index(
        proc::AbstractProcessDefinition,
        dir::ParticleDirection,
        species::AbstractParticleType,
        n::Int,
    )
    # find particle index of n-th particle of *this species and dir*
    species_index = 0
    for i in 1:n
        if particles(proc, dir)[i] == species
            species_index += 1
        end
    end

    return species_index
end

function _total_index(
        proc::AbstractProcessDefinition,
        dir::ParticleDirection,
        species::AbstractParticleType,
        n::Int,
    )
    # find particle index of all particles given n-th particle of dir and species (inverse of _species_index)
    total_index = 0
    species_count = 0
    for p in particles(proc, dir)
        total_index += 1
        if species == p
            species_count += 1
        end
        if species_count == n
            return if dir == Incoming()
                total_index
            else
                number_incoming_particles(proc) + total_index
            end
        end
    end

    throw("did not find $n-th $dir $species")
end

function _canonical_index(vp::VP) where {VP <: VirtualParticle}
    @assert vp.species != Photon "canonical index is only for (anti-)fermions"
    (left_ferms, right_ferms) = _canonical_fermion_indices(vp)

    # remove stuff
    for cycle in vp.open_cycles
        filter!(x -> x != cycle[1], left_ferms)
        filter!(x -> x != cycle[2], right_ferms)
    end

    left_minus_right = [setdiff(Set(left_ferms), Set(right_ferms))...]
    right_minus_left = [setdiff(Set(right_ferms), Set(left_ferms))...]

    if length(left_minus_right) == 1
        @assert isempty(right_minus_left)
        return (:left, left_minus_right[begin])
    elseif length(right_minus_left) == 1
        @assert isempty(left_minus_right)
        return (:right, right_minus_left[begin])
    else
        @assert false
    end
end

function _canonical_fermion_indices(vp::VP) where {VP <: VirtualParticle}
    left_ferms = Int[]
    right_ferms = Int[]
    for (contribs, dir) in Iterators.zip(
            (vp.in_particle_contributions, vp.out_particle_contributions),
            (Incoming(), Outgoing()),
        )
        c = 0
        for contrib in contribs
            c += 1
            if !contrib
                continue
            end
            (lr, index) = _get_canonical_index(vp.proc, dir, c)
            if lr == :left
                push!(left_ferms, index)
            elseif lr == :right
                push!(right_ferms, index)
            end
        end
    end

    return (left_ferms, right_ferms)
end
