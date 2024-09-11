#
# Feynman Diagram, tree-level, QED
#

"""
    FeynmanDiagram{N,E,U,T,M,FM} <: AbstractTreeLevelFeynmanDiagram

An implementation of [`AbstractTreeLevelFeynmanDiagram`](@ref), representing Feynman diagrams of tree-level perturbative QED.

The type parameters are:
- `N`: The total number of fermion lines in the diagram, i.e., `N := E + U + T`
- `E`: The total number of external electrons/positron pairs in the scattering process.
- `U`: The total number of external muon/antimuon pairs in the scattering process.
- `T`: The total number of external tauon/antitauon pairs in the scattering process.
- `M`: The total number of external photons in the scattering process.
- `FM`: A `FlatMatrix` type to efficiently store the diagram layout in serial memory.

!!! note
    `E`, `U`, and `T` are numbers of Fermion *lines* while `M` is the number of individual external photons. This means the number of external legs is `2(E + U + T) + M`.

!!! warning
    While `U` and `T` can be non-zero in this type and the rest of the code in this package is ready to deal with muons and tauons, no types for them exist yet in the QEDjl-project. Therefore, some functions will throw unimplemented errors when used with non-zero values for `U` and `T`.
"""
struct FeynmanDiagram{N,E,U,T,M,FM} <:
       AbstractTreeLevelFeynmanDiagram where {N,E,U,T,M,FM<:FlatMatrix}
    diagram_structure::FM

    electron_permutation::NTuple{E,Int}
    muon_permutation::NTuple{U,Int}
    tauon_permutation::NTuple{T,Int}

    function FeynmanDiagram(
        structure::Vector{Vector{Int}},
        elec_perm::Vector{Int},
        muon_perm::Vector{Int},
        tauon_perm::Vector{Int},
        ::Val{E},
        ::Val{U},
        ::Val{T},
        ::Val{M},
    ) where {E,U,T,M}
        @assert E == length(elec_perm)
        @assert U == length(muon_perm)
        @assert T == length(tauon_perm)
        N = E + U + T

        return new{N,E,U,T,M,FlatMatrix{Int64,2 * N - 2 + M,N}}(
            FlatMatrix(structure),
            NTuple{E,Int}(elec_perm),
            NTuple{U,Int}(muon_perm),
            NTuple{T,Int}(tauon_perm),
        )
    end
end

"""
    VirtualParticle{
        PROC<:AbstractProcessDefinition,
        PT<:AbstractParticleType,
        I,
        O,
    }

Representation of a virtual particle and the return type of the [`virtual_particles`](@ref) function.
The type parameters are:
- PROC: The process this particle is a process of.
- PT: The particle type of this virtual particle, e.g. [`QEDcore.Photon`](@extref) or [`QEDcore.Electron`](@extref).
- I: The number of incoming particles of the process, i.e. `I := length(incoming_paritcles(PROC()))`
- O: The number of outgoing particles of the process, i.e. `O := length(outgoing_paritcles(PROC()))`
"""
struct VirtualParticle{PROC<:AbstractProcessDefinition,PT<:AbstractParticleType,I,O}
    proc::PROC
    species::PT
    in_particle_contributions::NTuple{I,Bool}
    out_particle_contributions::NTuple{O,Bool}
end

function Base.show(io::IO, vp::VirtualParticle)
    pr = x -> x ? "1" : "0"
    return print(
        io,
        "$(vp.species): \t$(*(pr.(vp.in_particle_contributions)...)) | $(*(pr.(vp.out_particle_contributions)...))",
    )
end

"""
    process(::AbstractTreeLevelFeynmanDiagram)::QEDbase.AbstractProcessDefinition

Interface function that must be implemented for an instance of [`AbstractTreeLevelFeynmanDiagram`](@ref).

Return the specific [`QEDbase.AbstractProcessDefinition`](@extref) which the given diagram is for.
"""
function QEDbase.process(vp::VirtualParticle)
    return vp.proc
end

function QEDbase.particle_species(vp::VirtualParticle)
    return vp.species
end

_in_contributions(vp::VirtualParticle) = vp.in_particle_contributions
_out_contributions(vp::VirtualParticle) = vp.out_particle_contributions
_contributions(vp::VirtualParticle) = ((_in_contributions(vp), _out_contributions(vp)))

is_virtual(vp::VirtualParticle) = _number_contributions(vp) > 1
is_external(vp::VirtualParticle) = _number_contributions(vp) == 1

# "addition" of the bool tuples
# TODO: this should probably not overload and export a + operator for base types
function Base.:+(
    a::Tuple{NTuple{I,Bool},NTuple{O,Bool}}, b::Tuple{NTuple{I,Bool},NTuple{O,Bool}}
) where {I,O}
    # realistically, there should never be "colliding" 1s. if there are there is probably an error and this should be asserted
    #= for (i, j) in zip(a[1], b[1]) @assert !(i && j) end
    for (i, j) in zip(a[2], b[2]) @assert !(i && j) end =#

    return (ntuple(i -> a[1][i] || b[1][i], I), ntuple(i -> a[2][i] || b[2][i], O))
end

_invert(::Electron) = Positron()
_invert(::Positron) = Electron()
_invert(::Photon) = Photon()
function _invert(::AbstractParticleType)
    throw(InvalidInputError("unimplemented for this particle type"))
end

function _invert(virtual_particle::VirtualParticle)
    I = length(virtual_particle.in_particle_contributions)
    O = length(virtual_particle.out_particle_contributions)
    return VirtualParticle(
        virtual_particle.proc,
        _invert(virtual_particle.species),
        ntuple(x -> !virtual_particle.in_particle_contributions[x], I),
        ntuple(x -> !virtual_particle.out_particle_contributions[x], O),
    )
end

# normalize the representation
function normalize(virtual_particle::VirtualParticle)
    I = length(_in_contributions(virtual_particle))
    O = length(_out_contributions(virtual_particle))
    data = _contributions(virtual_particle)
    s = sum(data[1]) + sum(data[2])
    if s > (I + O) / 2
        return _invert(virtual_particle)
    elseif s == (I + O) / 2 && data[1][1] == false
        return _invert(virtual_particle)
    else
        return virtual_particle
    end
end

function _momentum_contribution(
    proc::AbstractProcessDefinition,
    dir::ParticleDirection,
    species::AbstractParticleType,
    index::Int,
)
    I = number_incoming_particles(proc)
    O = number_outgoing_particles(proc)

    # get index of n-th "dir species" particle in proc
    particles_seen = 0
    c = 0
    for p in particles(proc, dir)
        c += 1
        if p == species
            particles_seen += 1
        end
        if particles_seen == index
            return (
                ((is_incoming(dir) && x == c for x in 1:I)...,),
                ((is_outgoing(dir) && x == c for x in 1:O)...,),
            )
        end
    end

    @assert false "tried to get momentum contribution of $dir $species $index but it does not exist in $proc"
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
    if (n > 0 && n <= E)
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

function _external_particle(proc::AbstractProcessDefinition, n::Int)
    (dir, species, _) = _fermion_type(proc, n)
    if dir == Outgoing()
        species = _invert(species)
    end
    return VirtualParticle(proc, species, _momentum_contribution(proc, n)...)
end

function _number_contributions(vp::VirtualParticle)
    return sum(vp.in_particle_contributions) + sum(vp.out_particle_contributions)
end

function Base.isless(a::VirtualParticle, b::VirtualParticle)
    if _number_contributions(a) == _number_contributions(b)
        if a.in_particle_contributions == b.in_particle_contributions
            return a.out_particle_contributions < b.out_particle_contributions
        end
        return a.in_particle_contributions < b.in_particle_contributions
    end
    return _number_contributions(a) < _number_contributions(b)
end

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
    make_up(a::VirtualParticle, b::VirtualParticle, c::VirtualParticle)
    
For virtual particles `a`, `b`, and `c`, return true if `a` and `b`'s joint momentum contributions add up to `c`'s momentum contributions.
"""
function make_up(a::VirtualParticle, b::VirtualParticle, c::VirtualParticle)
    if particle_species(a) == Photon() && particle_species(b) == Photon()
        return false
    end
    # it should be unnecessary to check here that a and b can actually react. if a + b = c they must be able to if a, b and c all exist in the diagram.
    for (a_contrib, b_contrib, c_contrib) in
        Iterators.zip(Iterators.flatten.(_contributions.((a, b, c)))...)
        if c_contrib != a_contrib + b_contrib
            return false
        end
    end

    return true
end

"""
    are_total(a::VirtualParticle, b::VirtualParticle, c::VirtualParticle)

Return true if a, b and c combined contain all external particles exactly once.
"""
function are_total(a::VirtualParticle, b::VirtualParticle, c::VirtualParticle)
    for (a_contrib, b_contrib, c_contrib) in
        Iterators.zip(Iterators.flatten.(_contributions.((a, b, c)))...)
        if a_contrib + b_contrib + c_contrib != 1
            return false
        end
    end

    return true
end

function particle_pairs(particles::Vector)
    pairs = Dict{VirtualParticle,Vector{Tuple{VirtualParticle,VirtualParticle}}}()

    proc = process(first(particles))
    # make sure the "smallest" particles come first, i.e. those with few contributors
    all_particles = sort([_pseudo_virtual_particles(proc)..., particles...])

    # find pairs for every particle after the external ones (those can't have pairs)
    for p_i in
        (number_incoming_particles(proc) + number_outgoing_particles(proc) + 1):length(
        all_particles
    )
        p = all_particles[p_i]
        pairs[p] = Vector{Tuple{VirtualParticle,VirtualParticle}}()

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

function total_particle_triples(particles::Vector)
    # particle pairs making up the whole graph
    result_triples = Vector{Tuple{VirtualParticle,VirtualParticle,VirtualParticle}}()

    proc = process(first(particles))

    working_set = vcat(particles, _pseudo_virtual_particles(proc))

    photons = filter(p -> is_boson(particle_species(p)), working_set)

    # make electrons a set for fast deletion
    electrons = Set(
        filter(
            p -> is_fermion(particle_species(p)) && is_particle(particle_species(p)),
            working_set,
        ),
    )

    # make positrons a set for fast lookup
    positrons = Set(
        filter(
            p -> is_fermion(particle_species(p)) && is_anti_particle(particle_species(p)),
            working_set,
        ),
    )

    # no participant can have more than half the external particles, so every possible particle is contained here
    # every photon has exactly one electron and positron partner
    for ph in photons
        for e in electrons
            if !disjunct(ph, e)
                continue
            end

            # create the only partner the ph and e could have together, then look for it in the actual positrons
            expected_p = _invert(
                VirtualParticle(
                    proc, particle_species(e), (_contributions(ph) + _contributions(e))...
                ),
            )

            if expected_p in positrons
                @assert are_total(ph, e, expected_p)
                push!(result_triples, (ph, e, expected_p))
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

function virtual_particles(
    proc::AbstractProcessDefinition, diagram::FeynmanDiagram{N,E,U,T,M,FM}
) where {N,E,U,T,M,FM}
    fermion_lines = PriorityQueue{Int64,Int64}()

    # count number of internal photons in each fermion line and make a priority queue for fermion line => number of internal photons
    for i in 1:N
        count = 0
        for p in 1:length(diagram.diagram_structure, i)
            if diagram.diagram_structure[i, p] <= N
                # internal photon
                count += 1
            end
        end
        enqueue!(fermion_lines, i => count)
    end

    result = Vector()

    internal_photon_contributions = Dict()

    # 2: Loop: 
    # while there are incomplete fermion lines:
    #   take a fermion line where there is max=1 particle ∉ known_particles
    #   walk the fermion line, assign each virtual particle the momentum composition of the previous (or initial fermion if start) "+" the connected particle
    #   when/if the unknown particle is encountered, start walking from the other side
    #   when they meet at the unknown particle, assign the unknown particle Photon and left side - right side momentum contribution
    while !isempty(fermion_lines)
        current_line = dequeue!(fermion_lines)

        local unknown_photon_momentum = nothing
        # walk line from the *left* (either incoming electron or outgoing positron)
        (dir, species, _) = _fermion_type(proc, current_line)
        if dir == Outgoing()
            species = _invert(species)
        end

        cumulative_mom = _momentum_contribution(proc, current_line)

        for i in 1:length(diagram.diagram_structure, current_line)
            binding_particle = diagram.diagram_structure[current_line, i]
            if (binding_particle <= N) # binding_particle is an internal photon
                if haskey(internal_photon_contributions, binding_particle)   # if the binding particle is known
                    cumulative_mom += internal_photon_contributions[binding_particle]
                else # if the binding particle is unknown
                    # save so far momentum and break, add the right side momentum later
                    unknown_photon_momentum = cumulative_mom
                    break
                end
            else # binding_particle is an external photon
                cumulative_mom += _momentum_contribution(proc, binding_particle)
            end
            push!(result, VirtualParticle(proc, species, cumulative_mom...))
        end

        if isnothing(unknown_photon_momentum)
            # case where we're done (only one fermion line or last fermion line)
            # fermion_lines always has to be empty at this point, otherwise the tree wouldn't be connected
            @assert isempty(fermion_lines)
            continue
        end

        # find right side of the line
        right_line = diagram.electron_permutation[current_line]
        species = _invert(species)

        cumulative_mom = _momentum_contribution(proc, right_line)
        for i in length(diagram.diagram_structure, current_line):-1:1   # iterate from the right
            binding_particle = diagram.diagram_structure[current_line, i]
            if (binding_particle <= N) # binding_particle is an internal photon
                if haskey(internal_photon_contributions, binding_particle)   # if the binding particle is known, proceed as above
                    cumulative_mom += internal_photon_contributions[binding_particle]
                else # if the binding particle is unknown
                    # we have arrived at the "middle" of the line
                    # this line will be the unknown particle for the other lines
                    internal_photon_contributions[current_line] =
                        cumulative_mom + unknown_photon_momentum
                    # now we know that the fermion line that binding_particle binds to on the other end has one fewer unknown internal photons
                    fermion_lines[binding_particle] -= 1
                    # add the internal photon virtual particle
                    push!(
                        result,
                        VirtualParticle(
                            proc, Photon(), (cumulative_mom + unknown_photon_momentum)...
                        ),
                    )
                    break
                end
            else # binding_particle is an external photon
                cumulative_mom += _momentum_contribution(proc, binding_particle)
            end
            push!(result, VirtualParticle(proc, species, cumulative_mom...))
        end
    end

    return normalize.(result)[1:(end - 1)]
end

#
# Generate Feynman Diagrams
#

mutable struct ExternalPhotonIterator
    N::Int  # number of fermion lines
    M::Int  # number of external photons
    fermion_structures::Vector{Vector{Vector{Int}}}
    fermion_structure_index::Int
    photon_structures::Vector{Vector{Vector{Int}}}
    photon_structure_index::Int
end

function _feynman_structures(n::Int, m::Int)
    f = labelled_plane_trees(n)
    return ExternalPhotonIterator(n, m, f, 1, _external_photon_multiplicity(f[1], n, m), 1)
end

function Base.length(it::ExternalPhotonIterator)
    return factorial(it.M + 3 * it.N - 3, 2 * it.N - 1)
end

function Base.iterate(iter::ExternalPhotonIterator)
    return (iter.photon_structures[iter.photon_structure_index], nothing)
end

function Base.iterate(iter::ExternalPhotonIterator, n::Nothing)
    if iter.photon_structure_index == length(iter.photon_structures)
        if iter.fermion_structure_index == length(iter.fermion_structures)
            return nothing
        end

        iter.fermion_structure_index += 1
        iter.photon_structure_index = 1
        iter.photon_structures = _external_photon_multiplicity(
            iter.fermion_structures[iter.fermion_structure_index], iter.N, iter.M
        )
    else
        iter.photon_structure_index += 1
    end

    return (iter.photon_structures[iter.photon_structure_index], nothing)
end

function _external_photon_multiplicity(
    fermion_structure::Vector{Vector{Int}}, n::Int, m::Int
)
    if m == 0
        return [fermion_structure]
    end

    res = Vector{Vector{Vector{Int}}}()

    # go through lines
    for line_index in eachindex(fermion_structure)

        # go through indices start to end
        for index in 1:(length(fermion_structure[line_index]) + 1)
            # copy to prevent muting
            new_photon_structure = deepcopy(fermion_structure)
            # add new photon
            insert!(new_photon_structure[line_index], index, n + 1)
            # recurse
            append!(res, _external_photon_multiplicity(new_photon_structure, n + 1, m - 1))
        end
    end

    return res
end

mutable struct FeynmanDiagramIterator{E,U,T,M}
    e::Val{E}  # number of electron lines  (indices 1      - e)
    e_perms::Vector{Vector{Int}} # list of all the possible permutations of the electrons
    e_index::Int
    u::Val{U}  # number of muon lines      (indices e+1    - e+u)
    u_perms::Vector{Vector{Int}}
    u_index::Int
    t::Val{T}  # number of tauon lines     (indices e+u+1  - e+u+t)
    t_perms::Vector{Vector{Int}}
    t_index::Int
    m::Val{M}  # number of external photons
    photon_structure_iter::ExternalPhotonIterator
    photon_structure::Vector{Vector{Int}} # current structure that's being permuted
end

function Base.length(it::FeynmanDiagramIterator{E,U,T,M}) where {E,U,T,M}
    N = E + U + T
    return factorial(M + 3 * N - 3, 2 * N - 1) * factorial(E) * factorial(U) * factorial(T)
end

function Base.iterate(iter::FeynmanDiagramIterator{E,U,T,M}) where {E,U,T,M}
    N = E + U + T
    f = FeynmanDiagram(
        iter.photon_structure,
        iter.e_perms[iter.e_index],
        iter.u_perms[iter.u_index],
        iter.t_perms[iter.t_index],
        iter.e,
        iter.u,
        iter.t,
        iter.m,
    )
    return (f, nothing)
end

function Base.iterate(iter::FeynmanDiagramIterator{E,U,T,M}, ::Nothing) where {E,U,T,M}
    iter.t_index += 1

    if iter.t_index > length(iter.t_perms)
        iter.t_index = 1
        iter.u_index += 1
    end

    if iter.u_index > length(iter.u_perms)
        iter.u_index = 1
        iter.e_index += 1
    end

    if iter.e_index > length(iter.e_perms)
        iter.e_index = 1
        photon_iter_result = iterate(iter.photon_structure_iter, nothing)
        if isnothing(photon_iter_result)
            return nothing
        end
        (iter.photon_structure, _) = photon_iter_result
    end

    N = E + U + T
    f = FeynmanDiagram(
        iter.photon_structure,
        iter.e_perms[iter.e_index],
        iter.u_perms[iter.u_index],
        iter.t_perms[iter.t_index],
        iter.e,
        iter.u,
        iter.t,
        iter.m,
    )
    return (f, nothing)
end

"""
    feynman_diagrams(proc::AbstractProcessDefinition)

Return all tree-level Feynman diagrams that contribute to the given process, in perturbative QED.
"""
function feynman_diagrams(proc::AbstractProcessDefinition)
    return _feynman_diagrams(incoming_particles(proc), outgoing_particles(proc))
end

function _feynman_diagrams(in_particles::Tuple, out_particles::Tuple)
    count(x -> x isa Electron, in_particles) + count(x -> x isa Positron, out_particles) ==
    count(x -> x isa Positron, in_particles) + count(x -> x isa Electron, out_particles) ||
        throw(InvalidInputError("the given particles do not make a physical process"))

    # get the fermion line counts and external photon count
    e = count(x -> x isa Electron, in_particles) + count(x -> x isa Positron, out_particles)
    m = count(x -> x isa Photon, in_particles) + count(x -> x isa Photon, out_particles)
    # TODO: do this the same way as for e when muons and tauons are a part of QED.jl
    u = 0
    t = 0
    n = e + u + t

    # the numbers in the feynman diagram go as follows:
    # left electrons -> left muons -> left tauons -> left photons -> right photons -> right electrons -> right muons -> right tauons
    # a "left" fermion is simply an incoming fermion or outgoing antifermion of the type, while a "left" photon is an incoming photon, and the reverse for the right ones
    f_iter = _feynman_structures(e + u + t, m)
    e_perms = collect(permutations(Int[(n + m + 1):(n + m + e);]))
    u_perms = collect(permutations(Int[(n + m + e + 1):(n + m + e + u);]))
    t_perms = collect(permutations(Int[(n + m + e + u + 1):(n + m + e + u + t);]))
    first_photon_structure, _ = iterate(f_iter)

    return FeynmanDiagramIterator(
        Val(e),
        e_perms,
        1,
        Val(u),
        u_perms,
        1,
        Val(t),
        t_perms,
        1,
        Val(m),
        f_iter,
        first_photon_structure,
    )
end

# use a small LRU maxsize since these vectors could get large
@memoize LRU(maxsize=3) function virtual_particles(proc::AbstractProcessDefinition)
    # use a set for deduplication
    particles = Set{VirtualParticle}()
    for fd in feynman_diagrams(proc)
        push!(particles, virtual_particles(proc, fd)...)
    end

    # convert to vector
    return [particles...]
end
