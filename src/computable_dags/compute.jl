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

const e = sqrt(4π / 137)
_vertex() = -1im * e * gamma()

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

function compute( #=@inline=#
    ::ComputeTask_BaseState,
    input::BaseStateInput{PS,SPIN_POL},
) where {PS,SPIN_POL}
    species = particle_species(input.particle)
    if is_outgoing(input.particle)
        species = _invert(species)
    end
    return Propagated( # "propagated" because it goes directly into the next pair
        species,
        QEDbase.base_state(
            particle_species(input.particle),
            particle_direction(input.particle),
            momentum(input.particle),
            input.spin_pol,
        ),
        # bispinor, adjointbispinor, or lorentzvector
    )
end

struct PropagatorInput{VP_T<:VirtualParticle,PSP_T<:AbstractPhaseSpacePoint}
    vp::VP_T
    psp::Ref{PSP_T}
end

function compute( #=@inline=#
    ::ComputeTask_Propagator,
    input::PropagatorInput{VP_T,PSP_T},
) where {VP_T,PSP_T}
    vp_mom = zero(typeof(momentum(input.psp[], Incoming(), 1)))
    for i in eachindex(_in_contributions(input.vp))
        if _in_contributions(input.vp)[i]
            vp_mom += momentum(input.psp[], Incoming(), i)
        end
    end
    for o in eachindex(_out_contributions(input.vp))
        if (_out_contributions(input.vp))[o]
            vp_mom -= momentum(input.psp[], Outgoing(), o)
        end
    end

    vp_species = particle_species(input.vp)
    #println("$vp_species with $vp_mom")
    inner = QEDbase.propagator(vp_species, vp_mom)
    println("inner: $(inner[1, 1])")
    return inner
    # diracmatrix or scalar number
end

struct Unpropagated{PARTICLE_T<:AbstractParticleType,VALUE_T}
    particle::PARTICLE_T
    value::VALUE_T
end

function Base.:+(a::Unpropagated{P,V}, b::Unpropagated{P,V}) where {P,V}
    return Unpropagated(a.particle, a.value + b.value)
end

struct Propagated{PARTICLE_T<:AbstractParticleType,VALUE_T}
    particle::PARTICLE_T
    value::VALUE_T
end

# maybe add the γ matrix term here too?
function compute( #=@inline=#
    ::ComputeTask_Pair,
    electron::Propagated{Electron,V1},
    positron::Propagated{Positron,V2},
) where {V1,V2}
    return Unpropagated(Photon(), positron.value * _vertex() * electron.value)  # fermion - antifermion -> photon
end
function compute( #=@inline=#
    ::ComputeTask_Pair,
    positron::Propagated{Positron,V1},
    electron::Propagated{Electron,V2},
) where {V1,V2}
    return Unpropagated(Photon(), positron.value * _vertex() * electron.value)  # antifermion - fermion -> photon
end
function compute( #=@inline=#
    ::ComputeTask_Pair,
    photon::Propagated{Photon,V1},
    fermion::Propagated{F,V2},
) where {F<:FermionLike,V1,V2}
    return Unpropagated(fermion.particle, photon.value * _vertex() * fermion.value) # (anti-)fermion - photon -> (anti-)fermion
end
function compute( #=@inline=#
    ::ComputeTask_Pair,
    fermion::Propagated{F,V2},
    photon::Propagated{Photon,V1},
) where {F<:FermionLike,V1,V2}
    return Unpropagated(fermion.particle, photon.value * _vertex() * fermion.value) # photon - (anti-)fermion -> (anti-)fermion
end

function compute( #=@inline=#
    ::ComputeTask_PropagatePairs,
    left::PROP_V,
    right::Unpropagated{P,VAL},
) where {PROP_V,P<:AbstractParticleType,VAL}
    return Propagated(right.particle, left * right.value)
end
function compute( #=@inline=#
    ::ComputeTask_PropagatePairs,
    left::Unpropagated{P,VAL},
    right::PROP_V,
) where {PROP_V,P<:AbstractParticleType,VAL}
    return Propagated(left.particle, right * left.value)
end

function compute( #=@inline=#
    ::ComputeTask_Triple,
    photon::Propagated{Photon,V1},
    electron::Propagated{Electron,V2},
    positron::Propagated{Positron,V3},
) where {V1,V2,V3}
    return positron.value * _vertex() * photon.value * electron.value
end
function compute( #=@inline=#
    c::ComputeTask_Triple,
    photon::Propagated{Photon,V1},
    positron::Propagated{Positron,V2},
    electron::Propagated{Electron,V3},
) where {V1,V2,V3}
    return compute(c, photon, electron, positron)
end
function compute( #=@inline=#
    c::ComputeTask_Triple,
    f1::Propagated{F1,V1},
    f2::Propagated{F2,V2},
    photon::Propagated{Photon,V3},
) where {V1,V2,V3,F1<:FermionLike,F2<:FermionLike}
    return compute(c, photon, f1, f2)
end
function compute( #=@inline=#
    c::ComputeTask_Triple,
    f1::Propagated{F1,V1},
    photon::Propagated{Photon,V2},
    f2::Propagated{F2,V3},
) where {V1,V2,V3,F1<:FermionLike,F2<:FermionLike}
    return compute(c, photon, f1, f2)
end

# this compiles in a reasonable amount of time for up to about 1e4 parameters
# use a summation algorithm with more accuracy and/or parallelization
function compute(::ComputeTask_CollectPairs, args::Vararg{N,T}) where {N,T}
    println("summing $args")
    return sum(args)
end
function compute(::ComputeTask_CollectTriples, args::Vararg{N,T}) where {N,T}
    println("summing $args")
    return sum(args)
end
function compute(::ComputeTask_SpinPolCumulation, args::Vararg{N,T}) where {N,T} #=@inline=#
    sum = 0.0
    for arg in args
        sum += abs2(arg)
    end
    return sum
end
