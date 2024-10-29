struct ComputeTask_BaseState <: AbstractComputeTask end             # calculate the base state of an external particle
struct ComputeTask_Propagator <: AbstractComputeTask end            # calculate the propagator term of a virtual particle
struct ComputeTask_Pair <: AbstractComputeTask end                  # from a pair of virtual particle currents, calculate the product
struct ComputeTask_CollectPairs <: AbstractComputeTask              # for a list of virtual particle current pair products, sum
    children::Int
end
struct ComputeTask_PropagatePairs <: AbstractComputeTask end        # for the result of a CollectPairs compute task and a propagator, propagate the sum
struct ComputeTask_Triple <: AbstractComputeTask end                # from a triple of virtual particle currents, calculate the diagram result
struct ComputeTask_CollectTriples <: AbstractComputeTask            # sum over triples results and 
    children::Int
end
struct ComputeTask_SpinPolCumulation <: AbstractComputeTask         # abs2 sum over all spin/pol configs
    children::Int
end

# import so we don't have to repeat it all the time
import ComputableDAGs: compute, compute_effort, children

const e = sqrt(4π / 137.035999177)
const VERTEX = -1im * e * gamma()

compute_effort(::ComputeTask_BaseState) = 0
compute_effort(::ComputeTask_Propagator) = 0
compute_effort(::ComputeTask_Pair) = 0
compute_effort(::ComputeTask_CollectPairs) = 0
compute_effort(::ComputeTask_PropagatePairs) = 0
compute_effort(::ComputeTask_Triple) = 0
compute_effort(::ComputeTask_CollectTriples) = 0
compute_effort(::ComputeTask_SpinPolCumulation) = 0

children(::ComputeTask_BaseState) = 1
children(::ComputeTask_Propagator) = 1
children(::ComputeTask_Pair) = 2
children(t::ComputeTask_CollectPairs) = t.children
children(::ComputeTask_PropagatePairs) = 2
children(::ComputeTask_Triple) = 3
children(t::ComputeTask_CollectTriples) = t.children
children(t::ComputeTask_SpinPolCumulation) = t.children

struct BaseStateInput{PS_T<:AbstractParticleStateful,SPIN_POL_T<:AbstractSpinOrPolarization}
    particle::PS_T
    spin_pol::SPIN_POL_T
end

function compute(
    ::ComputeTask_BaseState, input::BaseStateInput{PS_T,SPIN_POL_T}
) where {PS_T<:AbstractParticleStateful,SPIN_POL_T<:AbstractSpinOrPolarization}
    species = particle_species(input.particle)
    if is_outgoing(input.particle)
        species = _invert(species)
    end
    state = QEDbase.base_state(
        particle_species(input.particle),
        particle_direction(input.particle),
        momentum(input.particle),
        input.spin_pol,
    )
    return Propagated( # "propagated" because it goes directly into the next pair
        species,
        state,
        # bispinor, adjointbispinor, or lorentzvector
    )
end

struct PropagatorInput{VP_T<:VirtualParticle,PSP_T<:AbstractPhaseSpacePoint}
    vp::VP_T
    psp::PSP_T
end

@inline _masked_sum(::Tuple{}, ::Tuple{}) = error("masked sum needs at least one argument")
@inline function _masked_sum(values::Tuple{T}, mask::Tuple{Bool}) where {T}
    return mask[1] ? values[1] : zero(T)
end
@inline function _masked_sum(
    values::Tuple{T,Vararg{T,N}}, mask::Tuple{Bool,Vararg{Bool,N}}
) where {N,T}
    return if mask[1]
        values[1] + _masked_sum(values[2:end], mask[2:end])
    else
        _masked_sum(values[2:end], mask[2:end])
    end
end

function _vp_momentum(
    vp::VirtualParticle{PROC,I,O}, psp::AbstractPhaseSpacePoint
) where {PROC,I,O}
    return _masked_sum(momenta(psp, Incoming()), _in_contributions(vp)) -
           _masked_sum(momenta(psp, Outgoing()), _out_contributions(vp))
end

function compute(
    ::ComputeTask_Propagator, input::PropagatorInput{VP_T,PSP_T}
) where {VP_T,PSP_T}
    vp_mom = _vp_momentum(input.vp, input.psp)
    vp_species = particle_species(input.vp)
    inner = QEDbase.propagator(vp_species, vp_mom)
    return inner
end

struct Unpropagated{PARTICLE_T<:AbstractParticleType,VALUE_T}
    particle::PARTICLE_T
    value::VALUE_T
end

@inline function Base.:+(a::Unpropagated{P,V}, b::Unpropagated{P,V}) where {P,V}
    return Unpropagated(a.particle, a.value + b.value)
end

struct Propagated{PARTICLE_T<:AbstractParticleType,VALUE_T}
    particle::PARTICLE_T
    value::VALUE_T
end

@inline function compute( # photon, electron
    ::ComputeTask_Pair,
    photon::Propagated{Photon},
    electron::Propagated{Electron},
)
    return Unpropagated(Electron(), photon.value * VERTEX * electron.value) # photon - electron -> electron
end
@inline function compute( # photon, positron
    ::ComputeTask_Pair,
    photon::Propagated{Photon},
    positron::Propagated{Positron},
)
    return Unpropagated(Positron(), positron.value * VERTEX * photon.value) # photon - positron -> positron
end
@inline function compute( # electron, positron
    ::ComputeTask_Pair,
    electron::Propagated{Electron},
    positron::Propagated{Positron},
)
    return Unpropagated(Photon(), positron.value * VERTEX * electron.value)  # electron - positron -> photon
end

@inline function compute(::ComputeTask_PropagatePairs, prop, photon::Unpropagated{Photon})
    return Propagated(Photon(), photon.value * prop)
end
@inline function compute(
    ::ComputeTask_PropagatePairs, prop, electron::Unpropagated{Electron}
)
    return Propagated(Electron(), prop * electron.value)
end
@inline function compute(
    ::ComputeTask_PropagatePairs, prop, positron::Unpropagated{Positron}
)
    return Propagated(Positron(), positron.value * prop)
end

@inline function compute(
    ::ComputeTask_Triple,
    photon::Propagated{Photon},
    electron::Propagated{Electron},
    positron::Propagated{Positron},
)
    return positron.value * (VERTEX * photon.value) * electron.value
end

# this compiles in a reasonable amount of time for up to about 1e4 parameters
# use a summation algorithm with more accuracy and/or parallelization
@inline function compute(::ComputeTask_CollectPairs, args::Vararg{N,T}) where {N,T}
    return sum(args)
end
@inline function compute(::ComputeTask_CollectTriples, args::Vararg{N,T}) where {N,T}
    return sum(args)
end
function compute(::ComputeTask_SpinPolCumulation, args::Vararg{N,T}) where {N,T}
    sum = 0.0
    for arg in args
        sum += abs2(arg)
    end
    return sum
end
