using Random
using QEDcore

# for a process definition, generate a random phase space point for the given process
# once QEDevents.jl supports random phase space point generation, that should be used here instead
function gen_process_input(rng::AbstractRNG, proc::AbstractProcessDefinition)
    in_momenta = ((rand_mom(rng, p) for p in incoming_particles(proc))...,)
    out_momenta = ((rand_mom(rng, p) for p in outgoing_particles(proc))...,)

    leftover = sum(in_momenta) - sum(out_momenta)

    local fermion_like_index
    for index in eachindex(out_momenta)
        if (incoming_particles(proc)[index] isa FermionLike)
            fermion_like_index = index
            break
        end
    end

    # this leads to the phase space point conserving momentum
    # it also makes the FermionLike off-shell
    out_momenta = ntuple(
        i -> i == fermion_like_index ? out_momenta[i] + leftover : out_momenta[i],
        length(out_momenta),
    )

    return PhaseSpacePoint(
        proc,
        PerturbativeQED(),
        PhasespaceDefinition(SphericalCoordinateSystem(), ElectronRestFrame()),
        in_momenta,
        out_momenta,
    )
end

# generate random on-shell photon
function rand_mom(rng::AbstractRNG, ::Photon)
    w = rand(rng)
    return SFourMomentum(w, 0.0, 0.0, w)
end

# generate random on-shell particle with mass
function rand_mom(rng::AbstractRNG, pt::AbstractParticleType)
    (x, y, z) = (rand(rng), rand(rng), rand(rng))
    E2 = x^2 + y^2 + z^2 + mass(pt)^2
    return SFourMomentum(sqrt(E2), x, y, z)
end
