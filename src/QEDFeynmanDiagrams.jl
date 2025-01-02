module QEDFeynmanDiagrams

using Reexport
@reexport using QEDbase
@reexport using QEDcore
@reexport using ComputableDAGs

using Combinatorics
using LRUCache
using Memoization
using DataStructures

export graph, number_of_diagrams

include("diagrams/virtual_particle.jl")
include("diagrams/vp_utils.jl")
include("diagrams/diagrams.jl")

include("computable_dags/compute.jl")
include("computable_dags/generation.jl")
include("computable_dags/fermion_sign.jl")
include("computable_dags/indexing.jl")
include("computable_dags/utils.jl")

end # module QEDFeynmanDiagrams
