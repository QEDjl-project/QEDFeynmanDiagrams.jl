_construction_string(::Incoming) = "Incoming()"
_construction_string(::Outgoing) = "Outgoing()"

_construction_string(::Electron) = "Electron()"
_construction_string(::Positron) = "Positron()"
_construction_string(::Photon) = "Photon()"

_construction_string(::PolX) = "PolX()"
_construction_string(::PolY) = "PolY()"
_construction_string(::SpinUp) = "SpinUp()"
_construction_string(::SpinDown) = "SpinDown()"

_species_str(::Photon) = "ph"
_species_str(::Electron) = "el"
_species_str(::Positron) = "po"

_spin_pol_str(::SpinUp) = "su"
_spin_pol_str(::SpinDown) = "sd"
_spin_pol_str(::PolX) = "px"
_spin_pol_str(::PolY) = "py"

_dir_str(::Incoming) = "inc"
_dir_str(::Outgoing) = "out"

# the possible spins or pols for generating base state tasks
_spin_pols(::AllSpin) = (SpinUp(), SpinDown())
_spin_pols(::SyncedSpin) = (SpinUp(), SpinDown())
_spin_pols(::SpinUp) = (SpinUp(),)
_spin_pols(::SpinDown) = (SpinDown(),)

_spin_pols(::AllPol) = (PolX(), PolY())
_spin_pols(::SyncedPol) = (PolX(), PolY())
_spin_pols(::PolX) = (PolX(),)
_spin_pols(::PolY) = (PolY(),)

_is_external(p::VirtualParticle) = _number_contributions(p) == 1

# return an index for the argument ordering on edges in the DAG for a given particle species, photon -> 1, electron -> 2, positron -> 3
_edge_index_from_species(::Photon) = 1
_edge_index_from_species(::Electron) = 2
_edge_index_from_species(::Positron) = 3
_edge_index_from_vp(vp::VirtualParticle) = _edge_index_from_species(particle_species(vp))

@inline _invert(::Electron) = Positron()
@inline _invert(::Positron) = Electron()
@inline _invert(::Photon) = Photon()

@inline _invert(t::Type) = typeof(_invert(t()))

function _invert(::AbstractParticleType)
    throw(InvalidInputError("unimplemented for this particle type"))
end

function _invert(virtual_particle::VirtualParticle)
    I = length(virtual_particle.in_particle_contributions)
    O = length(virtual_particle.out_particle_contributions)

    new_cycles = sort([(cycle[2], cycle[1]) for cycle in virtual_particle.open_cycles])

    return VirtualParticle(
        virtual_particle.proc,
        _invert(particle_species(virtual_particle)),
        ntuple(x -> !virtual_particle.in_particle_contributions[x], I),
        ntuple(x -> !virtual_particle.out_particle_contributions[x], O),
        new_cycles,
    )
end

Base.isless(::ParticleDirection, ::ParticleDirection) = false
Base.isless(::Incoming, ::Outgoing) = true
Base.isless(::UnknownDirection, ::Incoming) = true
Base.isless(::UnknownDirection, ::Outgoing) = true

function Base.isless(a::VirtualParticle, b::VirtualParticle)
    if _number_contributions(a) == _number_contributions(b)
        if a.in_particle_contributions == b.in_particle_contributions
            if a.out_particle_contributions == b.out_particle_contributions
                return a.open_cycles < b.open_cycles
            end
            return a.out_particle_contributions < b.out_particle_contributions
        end
        return a.in_particle_contributions < b.in_particle_contributions
    end
    return _number_contributions(a) < _number_contributions(b)
end
