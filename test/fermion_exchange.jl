# file for testing that fermion exchange negation is handled correctly

using Random
using QEDcore
using QEDbase.Mocks
using ComputableDAGs
using QEDFeynmanDiagrams

using Logging

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

include("utils.jl")
include("impl/bhabha.jl")

RNG = MersenneTwister(0)

@testset "Bhabha Scattering ep -> ep" begin
    g = graph(bhabha)
    # suppress type inference warnings; they don't matter here
    # they arise because the MockProcess isn't very type stable in its interface implementation
    # that leads to this being slow, but it's only mocking/testing, so that is okay
    f = with_logger(ConsoleLogger(Logging.Error)) do
        get_compute_function(g, bhabha, cpu_st(), @__MODULE__)
    end
    for i in 1:10
        input = gen_process_input(RNG, bhabha)

        @test isapprox(f(input), _ground_truth_bhabha(input))
    end
end
