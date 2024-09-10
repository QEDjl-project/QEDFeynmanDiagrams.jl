module QEDFeynmanDiagrams

using Reexport
@reexport using QEDbase
@reexport using QEDcore

using ComputableDAGs
using Combinatorics
using LRUCache
using Memoization
using DataStructures

export FeynmanDiagram, VirtualParticle
export feynman_diagrams
export external_particles, virtual_particles, process, generate_DAG

include("flat_matrix.jl")

include("diagrams/labelled_plane_trees.jl")
include("diagrams/interface.jl")
include("diagrams/diagrams.jl")

include("computable_dags/compute.jl")
include("computable_dags/generation.jl")

end # module QEDFeynmanDiagrams
