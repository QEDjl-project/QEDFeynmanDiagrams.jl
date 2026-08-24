Base.@irrational e (sqrt(4π / big(137.035999177)))

function VERTEX(::Type{F}) where {F <: AbstractFloat}
    return (-(one(complex(F))) * e * gamma(complex(F)))
end

struct BaseStateInput{PS_T <: AbstractParticleStateful, SPIN_POL_T <: AbstractSpinOrPolarization}
    particle::PS_T
    spin_pol::SPIN_POL_T

    function BaseStateInput(ps::PS_T, spinpol::SPIN_POL_T) where {PS_T, SPIN_POL_T}
        return new{PS_T, SPIN_POL_T}(ps, spinpol)
    end
end

struct PropagatorInput{VP_T <: VirtualParticle, PSP_T <: AbstractPhaseSpacePoint}
    vp::VP_T
    psp::PSP_T

    function PropagatorInput(vp::VP_T, psp::PSP_T) where {VP_T, PSP_T}
        return new{VP_T, PSP_T}(vp, psp)
    end
end

@inline _masked_sum(::Tuple{}, ::Tuple{}) = error("masked sum needs at least one argument")
@inline function _masked_sum(values::Tuple{T}, mask::Tuple{Bool}) where {T}
    return mask[1] ? values[1] : zero(T)
end
@inline function _masked_sum(
        values::Tuple{T, Vararg{T, N}}, mask::Tuple{Bool, Vararg{Bool, N}}
    ) where {N, T}
    return if mask[1]
        values[1] + _masked_sum(values[2:end], mask[2:end])
    else
        _masked_sum(values[2:end], mask[2:end])
    end
end

function _vp_momentum(
        vp::VirtualParticle{PROC, I, O}, psp::AbstractPhaseSpacePoint, ::Positron
    ) where {PROC, I, O}
    return -_masked_sum(momenta(psp, Incoming()), _in_contributions(vp)) +
        _masked_sum(momenta(psp, Outgoing()), _out_contributions(vp))
end

function _vp_momentum(
        vp::VirtualParticle{PROC, I, O}, psp::AbstractPhaseSpacePoint, ::AbstractParticleType
    ) where {PROC, I, O}
    return _masked_sum(momenta(psp, Incoming()), _in_contributions(vp)) -
        _masked_sum(momenta(psp, Outgoing()), _out_contributions(vp))
end

struct Unpropagated{PARTICLE_T <: AbstractParticleType, VALUE_T}
    particle::PARTICLE_T
    value::VALUE_T
end

@inline function Base.:+(a::Unpropagated{P, V}, b::Unpropagated{P, V}) where {P, V}
    return Unpropagated(a.particle, a.value + b.value)
end
@inline function Base.:*(z::Number, a::Unpropagated{P, V}) where {P, V}
    return Unpropagated(a.particle, z * a.value)
end
@inline function Base.:*(a::Unpropagated{P, V}, z::Number) where {P, V}
    return Unpropagated(a.particle, z * a.value)
end

struct Propagated{PARTICLE_T <: AbstractParticleType, VALUE_T}
    particle::PARTICLE_T
    value::VALUE_T
end

@compute_task ComputeTask_BaseState 0 function _base_state(input::BaseStateInput{PS_T, SPIN_POL_T}) where {PS_T <: AbstractParticleStateful, SPIN_POL_T <: AbstractSpinOrPolarization}
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

@compute_task ComputeTask_Propagator 0 function _propagator(input::PropagatorInput{VP_T, PSP_T}) where {VP_T, PSP_T}
    vp_species = particle_species(input.vp)
    vp_mom = _vp_momentum(input.vp, input.psp, vp_species)
    inner = QEDbase.propagator(vp_species, vp_mom)
    return inner
end

@compute_task ComputeTask_Pair 0 c_pair

function c_pair( # photon, electron
        photon::Propagated{Photon},
        electron::Propagated{Electron},
    )
    T = real(eltype(electron.value))
    return Unpropagated(Electron(), (photon.value * VERTEX(T)) * electron.value) # photon - electron -> electron
end
function c_pair( # photon, positron
        photon::Propagated{Photon},
        positron::Propagated{Positron},
    )
    T = real(eltype(positron.value))
    return Unpropagated(Positron(), positron.value * (VERTEX(T) * photon.value)) # photon - positron -> positron
end
function c_pair( # electron, positron
        electron::Propagated{Electron},
        positron::Propagated{Positron},
    )
    T = real(eltype(electron.value))
    return Unpropagated(Photon(), positron.value * VERTEX(T) * electron.value)  # electron - positron -> photon
end

@compute_task ComputeTask_PairNegated 0 function _pair_negated(v1::Propagated{P1}, v2::Propagated{P2}) where {P1, P2}
    T = real(eltype(v1.value))
    return -one(T) * c_pair(v1, v2)
end

@compute_task ComputeTask_PropagatePairs 0 (
    c_prop_pairs(prop, photon::Unpropagated{Photon}) = Propagated(Photon(), photon.value * prop);
    c_prop_pairs(prop, electron::Unpropagated{Electron}) = Propagated(Electron(), prop * electron.value);
    c_prop_pairs(prop, positron::Unpropagated{Positron}) = Propagated(Positron(), positron.value * prop)
)

@compute_task ComputeTask_Triple 0 function _triple(
        photon::Propagated{Photon},
        electron::Propagated{Electron},
        positron::Propagated{Positron},
    )
    T = real(eltype(photon.value))
    return positron.value * (VERTEX(T) * photon.value) * electron.value
end

@compute_task ComputeTask_TripleNegated 0 function _triple_negated(
        photon::Propagated{Photon},
        electron::Propagated{Electron},
        positron::Propagated{Positron},
    )
    T = real(eltype(photon.value))
    return -one(T) * (positron.value * (VERTEX(T) * photon.value) * electron.value)
end

# this compiles in a reasonable amount of time for up to about 1e4 parameters
# TODO: use a summation algorithm with more accuracy and/or parallelization
@compute_task ComputeTask_CollectPairs 0 function _sum_pairs(args::Vararg)
    sum(args)
end

@compute_task ComputeTask_CollectTriples 0 function _sum_triples(args::Vararg)
    sum(args)
end

@compute_task ComputeTask_SpinPolCumulation 0 function _sum_spin_pol(args::Vararg)
    sum(abs2, args)
end

# for differential probability and cross-sections overloads
@compute_task ComputeTask_UnsafeDiffProb 0 function _diff_prob(mat_el_sqsum::T, psp) where {T}
    normalization = QEDbase._averaging_norm(T, psp.proc)
    ps_fac = QEDbase._phase_space_factor(psp)
    return normalization * mat_el_sqsum * ps_fac
end

@compute_task ComputeTask_UnsafeDiffCS 0 function _diff_cs(diff_prob::T, psp) where {T}
    return 1 / (4 * QEDbase._incident_flux(psp)) * diff_prob
end
