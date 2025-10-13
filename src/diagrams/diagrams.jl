"""
    number_of_diagrams(proc::AbstractProcessDefinition)

For a given [`QEDbase.AbstractProcessDefinition`](@extref), returns the number of valid
Feynman diagrams at tree-level. This is equivalent to

\$\\frac{(M + 3N - 3)!}{(2N - 1)!} * E! * U! * T!\$\\
, where\\
\$M \\dots\$ number of external photons,\\
\$E \\dots\$ number of electron lines,\\
\$U \\dots\$ number of muon lines,\\
\$T \\dots\$ number of tauon lines, and\\
\$N = E + U + T\$.\\

An electron/muon/tauon "line" is a pair of incoming fermion and outgoing anti-fermion or vice versa.
"""
function number_of_diagrams(proc::AbstractProcessDefinition)
    M =
        number_particles(proc, Incoming(), Photon()) +
        number_particles(proc, Outgoing(), Photon())

    E =
        number_particles(proc, Incoming(), Electron()) +
        number_particles(proc, Outgoing(), Positron())
    anti_E =
        number_particles(proc, Incoming(), Positron()) +
        number_particles(proc, Outgoing(), Electron())

    if E != anti_E || (E == anti_E == 0)
        return 0
    end

    # TODO: add muons/tauons
    U = 0
    T = 0

    N = E + U + T
    return factorial(M + 3 * N - 3, 2 * N - 1) * factorial(E) * factorial(U) * factorial(T)
end

# "addition" of the bool tuples
# TODO: this should probably not overload and export a + operator for base types
function Base.:+(
        a::Tuple{NTuple{I, Bool}, NTuple{O, Bool}}, b::Tuple{NTuple{I, Bool}, NTuple{O, Bool}}
    ) where {I, O}
    return (ntuple(i -> a[1][i] != b[1][i], I), ntuple(i -> a[2][i] != b[2][i], O))
end

@inline function _momentum_contribution_helper(
        proc::AbstractProcessDefinition,
        parts::Tuple{},
        dir::ParticleDirection,
        species::AbstractParticleType,
        index::Int,
        c::Int,
    )
    throw(
        "tried to get momentum contribution of $species $index but it does not exist in $proc",
    )
end
@inline function _momentum_contribution_helper(
        proc::AbstractProcessDefinition,
        parts::Tuple{SPECIES1, Vararg},
        species::SPECIES2,
        dir::DIR,
        index::Int,  # index of particle to find
        c::Int,      # count of seen particles
    ) where {SPECIES1, SPECIES2, DIR}
    return _momentum_contribution_helper(proc, parts[2:end], species, dir, index, c + 1)
end
@inline function _momentum_contribution_helper(
        proc::AbstractProcessDefinition,
        parts::Tuple{SPECIES, Vararg},
        species::SPECIES,
        dir::DIR,
        index::Int, # index of particle to find
        c::Int,      # count of seen particles
    ) where {DIR, SPECIES}
    # equal species, check index, then call next
    if index == 0
        return (
            ntuple(x -> (is_incoming(dir) && x == c), number_incoming_particles(proc)),
            ntuple(x -> (is_outgoing(dir) && x == c), number_outgoing_particles(proc)),
        )
    end
    return _momentum_contribution_helper(proc, parts[2:end], species, dir, index - 1, c + 1)
end

function _momentum_contribution(
        proc::AbstractProcessDefinition,
        dir::ParticleDirection,
        species::AbstractParticleType,
        index::Int,
    )
    return _momentum_contribution_helper(
        proc, particles(proc, dir), species, dir, index - 1, 1
    )
end

function _fermion_type(proc::AbstractProcessDefinition, n::Int)
    E =
        number_particles(proc, Incoming(), Electron()) +
        number_particles(proc, Outgoing(), Positron())
    U = 0 # TODO add muons
    T = 0 # TODO add tauons
    M =
        number_particles(proc, Incoming(), Photon()) +
        number_particles(proc, Outgoing(), Photon())
    N = E + U + T

    # from the fermion index, get (Direction, Species, n) tuple, where n means it's the nth particle of that dir and species
    return if (n > 0 && n <= E)
        electron_n = n
        if electron_n > number_particles(proc, Incoming(), Electron())
            return (
                Outgoing(),
                Positron(),
                electron_n - number_particles(proc, Incoming(), Electron()),
            )
        else
            return (Incoming(), Electron(), electron_n)
        end
    elseif (n > E && n <= E + U)
        # left muon n - E
        muon_n = n - E
        throw(InvalidInputError("unimplemented for muons"))
    elseif (n > E + U && n <= E + U + T)
        # left tauon n - E - U
        tauon_n = n - E - U
        throw(InvalidInputError("unimplemented for tauons"))
    elseif (n > N && n <= N + M)
        # photon
        photon_n = n - N
        if photon_n > number_particles(proc, Incoming(), Photon())
            return (
                Outgoing(),
                Photon(),
                photon_n - number_particles(proc, Incoming(), Photon()),
            )
        else
            return (Incoming(), Photon(), photon_n)
        end
    elseif (n > N + M && n <= N + M + E)
        # right electron
        electron_n = n - N - M
        if electron_n > number_particles(proc, Outgoing(), Electron())
            # incoming positron
            return (
                Incoming(),
                Positron(),
                electron_n - number_particles(proc, Outgoing(), Electron()),
            )
        else
            # outgoing electron
            return (Outgoing(), Electron(), electron_n)
        end
    elseif (n > N + M + E && n <= N + M + E + U)
        # right muon
        muon_n = n - N - M - E
        throw(InvalidInputError("unimplemented for muons"))
    elseif (n > N + M + E + U && n <= N + M + E + U + T)
        # right tauon
        tauon_n = n - N - M - E - U
        throw(InvalidInputError("unimplemented for tauons"))
    else
        # error
        throw(InvalidInputError("invalid index given"))
    end
end

@inline function _momentum_contribution(proc::AbstractProcessDefinition, n::Int)
    return _momentum_contribution(proc, _fermion_type(proc, n)...)
end

function _external_particle(proc::PROC, n::Int) where {PROC <: AbstractProcessDefinition}
    I = number_incoming_particles(proc)
    O = number_outgoing_particles(proc)
    SPECIFIC_VP = VirtualParticle{PROC, NTuple{I, Bool}, NTuple{O, Bool}}

    (dir, species, _) = _fermion_type(proc, n)
    if dir == Outgoing()
        species = _invert(species)
    end
    return SPECIFIC_VP(proc, species, _momentum_contribution(proc, n)...)
end

function _number_contributions(vp::VirtualParticle)
    return sum(vp.in_particle_contributions) + sum(vp.out_particle_contributions)
end

"""
    particle_pairs(particles::Vector{VirtualParticle})

From a vector of particles (e.g., generated from [`virtual_particles`](@ref)), generate
a `Dict` which maps from a [`VirtualParticle`](@ref) to a vector of `Tuple`s (pairs) of
[`VirtualParticle`](@ref)s. The two virtual particles of each tuple [`make_up`](@ref) the
key `VirtualParticle`.

The result is used in the [`graph`](@ref) generation.

See also: [`total_particle_triples`](@ref)
"""
function particle_pairs(
        particles::Vector{SPECIFIC_VP}
    ) where {PROC, I, O, SPECIFIC_VP <: VirtualParticle{PROC, I, O}}
    pairs = Dict{SPECIFIC_VP, Vector{Tuple{SPECIFIC_VP, SPECIFIC_VP}}}()

    proc = process(first(particles))
    # make sure the "smallest" particles come first, i.e. those with few contributors
    all_particles::Vector{SPECIFIC_VP} = _pseudo_virtual_particles(proc)
    append!(all_particles, particles)
    sort!(all_particles)

    # find pairs for every particle after the external ones (those can't have pairs)
    for p_i in
        (number_incoming_particles(proc) + number_outgoing_particles(proc) + 1):length(
            all_particles
        )
        p = all_particles[p_i]
        pairs[p] = Vector{Tuple{SPECIFIC_VP, SPECIFIC_VP}}()

        # only need to consider external particles and virtual particles that come before p_i
        for p_a_i in 1:(p_i - 2)
            # and only partners between a and p_i
            for p_b_i in (p_a_i + 1):(p_i - 1)
                p_a = all_particles[p_a_i]
                p_b = all_particles[p_b_i]

                if make_up(p_a, p_b, p)
                    push!(pairs[p], (p_a, p_b))
                end
            end
        end
    end

    return pairs
end

"""
    total_particle_triples(particles::Vector{VirtualParticle})

Similar to [`particle_pairs`](@ref), this generates a `Vector` of `Tuples` (triples). Each tuple
contains three particles, a [`QEDcore.Photon`](@extref), a [`QEDcore.Fermion`](@extref), and a
[`QEDcore.AntiFermion`](@extref). These three particles [`are_total`](@ref).
"""
function total_particle_triples(
        particles::Vector{VirtualParticle{PROC, I, O}}
    ) where {PROC, I, O}
    SPECIFIC_VP = VirtualParticle{PROC, I, O}
    # particle pairs making up the whole graph
    result_triples = Vector{Tuple{SPECIFIC_VP, SPECIFIC_VP, SPECIFIC_VP}}()

    proc = process(first(particles))

    working_set = vcat(particles, _pseudo_virtual_particles(proc))

    photons = filter(p -> particle_species(p) == Photon(), working_set)

    # make electrons a set for fast deletion
    electrons = Set(filter(p -> particle_species(p) == Electron(), working_set))

    # make positrons a set for fast lookup
    positrons = Set(filter(p -> particle_species(p) == Positron(), working_set))

    # no participant can have more than half the external particles, so every possible particle is contained here
    # every photon has exactly one electron and positron partner
    for ph in photons
        for e in electrons
            if !disjunct(ph, e)
                continue
            end

            for p in positrons
                if are_total(ph, e, p)
                    push!(result_triples, (ph, e, p))
                end
            end
        end
    end

    return result_triples
end

"""
    _pseudo_virtual_particles

Return a vector of [`VirtualParticle`](@ref) for each external particle. These are not actually virtual particles, but can be helpful as entry points.
"""
function _pseudo_virtual_particles(proc::AbstractProcessDefinition)
    return sort(
        _external_particle.(
            proc, [1:(number_incoming_particles(proc) + number_outgoing_particles(proc));]
        ),
    )
end

"""
    _count_particles(particles::Tuple{SPECIES...}, ::Tuple{Bool...}, species::AbstractParticleType)

Return the number of particles of the given species in the first tuple where at the same index the second tuple contains a `true`.

The tuples need to match in length or an error is thrown.
"""
@inline _count_particles(::Tuple{}, ::Tuple{}, species) = 0
@inline function _count_particles(
        parts::Tuple{SPECIES, Vararg}, bools::Tuple{Bool, Vararg}, species::SPECIES
    ) where {SPECIES}
    return (bools[1] ? 1 : 0) + _count_particles(parts[2:end], bools[2:end], species)
end
@inline function _count_particles(
        parts::Tuple{SPECIES1, Vararg}, bools::Tuple{Bool, Vararg}, species::SPECIES2
    ) where {SPECIES1, SPECIES2}
    return 0 + _count_particles(parts[2:end], bools[2:end], species)
end

function reduce_cycles(vec::Vector{OPEN_FERMION_CYCLE_T})
    # CAUTION: chat-gpt generated, but tested
    while true
        # Flag to check if changes occur
        changed = false
        new_vec = OPEN_FERMION_CYCLE_T[]
        skip_indices = Set{Int}()

        for i in 1:length(vec)
            if i in skip_indices
                continue
            end

            fused = false
            for j in (i + 1):length(vec)
                if j in skip_indices
                    continue
                end

                a, b = vec[i]
                c, d = vec[j]

                # Fuse if b == c
                if b == c
                    push!(new_vec, (a, d))
                    push!(skip_indices, j)
                    changed = true
                    fused = true
                    break
                    # Fuse if d == a
                elseif d == a
                    push!(new_vec, (c, b))
                    push!(skip_indices, j)
                    changed = true
                    fused = true
                    break
                end
            end

            # Add the original tuple if it wasn't fused and numbers aren't equal
            if !fused && vec[i][1] != vec[i][2]
                push!(new_vec, vec[i])
            end
        end

        # Update the vector
        vec = new_vec

        # Break if no changes
        if !changed
            break
        end
    end

    return vec
end

function _count_closed_cycles(vec::Vector{OPEN_FERMION_CYCLE_T})
    closed_count = 0
    # Adapted from reduce_cycles
    while true
        # Flag to check if changes occur
        changed = false
        new_vec = OPEN_FERMION_CYCLE_T[]
        skip_indices = Set{Int}()

        for i in 1:length(vec)
            if i in skip_indices
                continue
            end

            fused = false
            for j in (i + 1):length(vec)
                if j in skip_indices
                    continue
                end

                a, b = vec[i]
                c, d = vec[j]

                # Fuse if b == c
                if b == c
                    push!(new_vec, (a, d))
                    push!(skip_indices, j)
                    changed = true
                    fused = true
                    break
                    # Fuse if d == a
                elseif d == a
                    push!(new_vec, (c, b))
                    push!(skip_indices, j)
                    changed = true
                    fused = true
                    break
                end
            end

            # Add the original tuple if it wasn't fused and numbers aren't equal
            if !fused
                if vec[i][1] == vec[i][2]
                    closed_count += 1   # only here do we actually close a cycle
                else
                    push!(new_vec, vec[i])
                end
            end
        end

        # Update the vector
        vec = new_vec

        # Break if no changes
        if !changed
            break
        end
    end

    return closed_count
end

function _find_cycles(left_ferms::Vector{Int}, right_ferms::Vector{Int})
    all_cycles = OPEN_FERMION_CYCLE_T[]
    for (l, r) in Iterators.zip(left_ferms, right_ferms)
        push!(all_cycles, (l, r))
    end

    return reduce_cycles(all_cycles)
end

# TODO: @memoize this?
function _open_cycle_helper(left_ferms::Vector{Int}, right_ferms::Vector{Int})
    result = Set{Vector{OPEN_FERMION_CYCLE_T}}()

    n = min(length(left_ferms), length(right_ferms))

    if (length(left_ferms) == n)
        for right_ferm_perm in permutations(right_ferms)
            cycles = _find_cycles(left_ferms, right_ferm_perm[1:n])
            push!(result, cycles)
        end
    else
        for left_ferm_perm in permutations(left_ferms)
            cycles = _find_cycles(left_ferm_perm[1:n], right_ferms)
            push!(result, cycles)
        end
    end

    return sort([result...])
end

"""
    gen_specific_vp_with_open_cycles

Return a `Vector` of all possible virtual particles with this configuration, i.e., with all possible distinct open cycles.
"""
function gen_specific_vp_with_open_cycles(
        proc::PROC, species::SPECIES, in_contribs::NTuple{I, Bool}, out_contribs::NTuple{O, Bool}
    ) where {PROC <: AbstractProcessDefinition, SPECIES <: AbstractParticleType, I, O}
    # get canonical indices of all participating fermions
    left_ferms = Int[]
    right_ferms = Int[]
    for (contribs, dir) in
        Iterators.zip((in_contribs, out_contribs), (Incoming(), Outgoing()))
        c = 0
        for contrib in contribs
            c += 1
            if !contrib
                continue
            end
            (lr, index) = _get_canonical_index(proc, dir, c)
            if lr == :left
                push!(left_ferms, index)
            elseif lr == :right
                push!(right_ferms, index)
            end
        end
    end

    open_cycles = _open_cycle_helper(left_ferms, right_ferms)
    return [
        VirtualParticle(proc, species, in_contribs, out_contribs, open_cycle) for
            open_cycle in open_cycles
    ]
end

"""
    virtual_particles(proc::AbstractProcessDefinition)

For a given [`QEDbase.AbstractProcessDefinition`](@extref), generate all virtual particles ([`VirtualParticle`](@ref)) that occur
in some valid diagram. For more information see the virtual particle docs.
"""
@memoize LRU(maxsize = 3) function virtual_particles(
        proc::PROC
    ) where {PROC <: AbstractProcessDefinition}
    I = number_incoming_particles(proc)
    O = number_outgoing_particles(proc)

    total_electrons =
        number_particles(proc, Incoming(), Electron()) +
        number_particles(proc, Outgoing(), Positron())

    SPECIFIC_VP = VirtualParticle{PROC, NTuple{I, Bool}, NTuple{O, Bool}}
    # use a set for deduplication
    particles = SPECIFIC_VP[]

    in_p = incoming_particles(proc)
    out_p = outgoing_particles(proc)

    for (in_contribs, out_contribs) in Iterators.product(
            Iterators.product(ntuple(_ -> (false, true), I)...),
            Iterators.product(ntuple(_ -> (false, true), O)...),
        )
        # check whether the contributions make a valid particle
        electrons =
            _count_particles(in_p, in_contribs, Electron()) +
            _count_particles(out_p, out_contribs, Positron())
        positrons =
            _count_particles(in_p, in_contribs, Positron()) +
            _count_particles(out_p, out_contribs, Electron())
        photons =
            _count_particles(in_p, in_contribs, Photon()) +
            _count_particles(out_p, out_contribs, Photon())

        # sort out invalid combinations
        if electrons + positrons + photons <= 1
            continue
        elseif electrons + positrons + photons > (I + O) / 2
            continue
        elseif electrons + positrons + photons == (I + O) / 2 && in_contribs[1] == false
            continue
        end

        # infer the species type
        local species::Type
        if electrons - positrons == 1
            # regardless of number of photons, one electron "leftover" means the result is a electron
            species = Electron
        elseif positrons - electrons == 1
            # regardless of number of photons, one positron "leftover" means the result is a positron
            species = Positron
        elseif electrons == positrons
            if electrons == 0 && photons > 1
                # multiple photons cannot interact without an electron or positron
                continue
            end
            if electrons == total_electrons
                # cannot "use up" all electrons and have photons leftover
                continue
            end

            species = Photon
        else
            continue
        end

        vps = gen_specific_vp_with_open_cycles(proc, species(), in_contribs, out_contribs)
        for vp in vps
            push!(particles, vp)
        end
    end
    return particles
end
