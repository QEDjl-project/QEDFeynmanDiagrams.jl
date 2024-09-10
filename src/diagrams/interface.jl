"""
    AbstractTreeLevelFeynmanDiagram

Abstract base type for FeynmanDiagrams. Must implement the functions

```julia
process(::AbstractTreeLevelFeynmanDiagram)::QEDbase.AbstractProcessDefinition
virtual_particles(::AbstractTreeLevelFeynmanDiagram)::NTuple{N, Tuple{QEDbase.AbstractParticleType, BitArray}}
```

By using the `QEDbase.AbstractProcessDefinition` interface, the function [`external_particles`](@ref) is automatically provided.

For more information on what the interface functions should do, see their documentation: [`process`](@ref), [`virtual_particles`](@ref)
"""
abstract type AbstractTreeLevelFeynmanDiagram end

"""
    virtual_particles(::QEDbase.AbstractProcessDefinition, ::AbstractTreeLevelFeynmanDiagram)::Vector{VirtualParticle}

Interface function that must be implemented for an instance of [`AbstractTreeLevelFeynmanDiagram`](@ref).

Return an `NTuple` with N elements, where N is the number of virtual particles in this diagram. For tree-level Feynman diagrams, \$N = k - 3\$, where \$k\$ is the number of external particles.
The elements of the `NTuple` are themselves `Tuple`s, containing for each virtual particle its `QEDbase.AbstractParticleType` and an `NTuple{, Bool}` indicating
with a `1` that an incoming external particle's momentum contributes to the virtual particle's momentum, and a `0` otherwise. The second `NTuple{, Bool}` does the same for the outgoing external
particles, which contribute their momentum negatively.
From this definition follows that a particles' `Boolean NTuple`ss are equivalent to their inverse, i.e., `BitArray`s where every bit is negated.

Example: Consider the Compton scattering process \$e^- + \\gamma \\to e^- + \\gamma\$ with the diagram where the incoming electron interacts with the incoming photon first.
For this diagram there is exactly one virtual particle, which is an electron. This electron's momentum can be represented as the sum of the two incoming particles' momenta, or 
that of the two outgoing particles. In the second possible diagram, where the incoming electron interacts with the outgoing photon first, the virtual particle is still an electron
but its momentum is the sum of the momenta of the incoming electron and the outgoing photon, or, equivalently, the outgoing electron and the incoming photon.


    virtual_particles(::AbstractProcessDefinition)::Vector{VirtualParticle}

Function that returns all unique virtual particles of the given process.

!!! note
    This function is usually costly to compute and used across multiple functions. Therefore, it caches its results using `Memoization.jl`.
"""
function virtual_particles end

"""
    external_particles(diagram::AbstractTreeLevelFeynmanDiagram)

Return a tuple of the incoming and outgoing particles (`QEDbase.AbstractParticleType`) of the diagram.
"""
function external_particles(diagram::AbstractTreeLevelFeynmanDiagram)
    return (incoming_particles(process(diagram)), outgoing_particles(process(diagram)))
end
