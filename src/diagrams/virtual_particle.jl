const OPEN_FERMION_CYCLE_T = Tuple{Int64,Int64}

"""
    VirtualParticle{
        PROC<:AbstractProcessDefinition,
        NTuple{I,Bool},
        NTuple{O,Bool},
    }

Representation of a virtual particle and the return type of the [`virtual_particles`](@ref) function.
The type parameters are:
- PROC: The process this particle is a process of.
- PT: The particle type of this virtual particle, e.g. [`QEDcore.Photon`](@extref) or [`QEDcore.Electron`](@extref).
- I: NTuple of Bools with the incoming momentum contributions
- O: NTuple of Bools with the outgoing momentum contributions
"""
struct VirtualParticle{PROC<:AbstractProcessDefinition,IT<:NTuple,OT<:NTuple}
    proc::PROC
    species::Type
    in_particle_contributions::IT
    out_particle_contributions::OT

    # open cycles in the context of fermion permutations
    # for n fermion lines in a process, there can be between 1 (like 1-2, 2-3, 3-1) and n (like 1-1, 2-2, 3-3) cycles
    # where the left number represents the canonical fermion index and the right number the canonical antifermion index
    open_cycles::Vector{OPEN_FERMION_CYCLE_T}

    function VirtualParticle(
        proc::PROC, species::PT, in_contrib::I, out_contrib::O
    ) where {PROC,PT,I,O}
        return new{PROC,I,O}(
            proc, typeof(species), in_contrib, out_contrib, OPEN_FERMION_CYCLE_T[]
        )
    end
    function VirtualParticle{PROC,I,O}(
        proc::PROC, species::PT, in_contrib::I, out_contrib::O
    ) where {PROC,PT,I,O}
        return new{PROC,I,O}(
            proc, typeof(species), in_contrib, out_contrib, OPEN_FERMION_CYCLE_T[]
        )
    end
    function VirtualParticle(
        proc::PROC,
        species::PT,
        in_contrib::I,
        out_contrib::O,
        open_cycles::Vector{OPEN_FERMION_CYCLE_T},
    ) where {PROC,PT,I,O}
        return new{PROC,I,O}(proc, typeof(species), in_contrib, out_contrib, open_cycles)
    end
end

function Base.hash(vp::VP, h::UInt) where {VP<:VirtualParticle}
    h = hash(VP, h)
    h = hash(vp.proc, h)
    h = hash(vp.species, h)
    h = hash(vp.in_particle_contributions, h)
    h = hash(vp.out_particle_contributions, h)
    h = hash(vp.open_cycles, h)
    return h
end

function Base.isequal(vp1::VP, vp2::VP) where {VP<:VirtualParticle}
    return vp1.species == vp2.species &&
           vp1.in_particle_contributions == vp2.in_particle_contributions &&
           vp1.out_particle_contributions == vp2.out_particle_contributions &&
           vp1.open_cycles == vp2.open_cycles
end

function Base.show(io::IO, vp::VirtualParticle)
    pr = x -> x ? "1" : "0"
    return print(
        io,
        "$(string(particle_species(vp))[1:3]): $(*(pr.(vp.in_particle_contributions)...)) | $(*(pr.(vp.out_particle_contributions)...)) | $(isempty(vp.open_cycles) ? "[      ]" : "$(vp.open_cycles)")",
    )
end

@inline function QEDbase.process(vp::VirtualParticle)
    return vp.proc
end

@inline function QEDbase.particle_species(vp::VirtualParticle)
    return (vp.species)()
end

@inline function _in_contributions(vp::VirtualParticle{PROC,I,O})::I where {PROC,I,O}
    return vp.in_particle_contributions
end
@inline function _out_contributions(vp::VirtualParticle{PROC,I,O})::O where {PROC,I,O}
    return vp.out_particle_contributions
end
@inline function _contributions(vp::VirtualParticle{PROC,I,O})::Tuple{I,O} where {PROC,I,O}
    return ((_in_contributions(vp), _out_contributions(vp)))
end

@inline function is_virtual(vp::VirtualParticle)
    return _number_contributions(vp) > 1
end
@inline function is_external(vp::VirtualParticle)
    return _number_contributions(vp) == 1
end
