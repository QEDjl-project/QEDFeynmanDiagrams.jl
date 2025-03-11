# file for testing that fermion exchange negation is handled correctly

using Random
using QEDcore
using QEDprocesses
using ComputableDAGs
using QEDFeynmanDiagrams

using Logging

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

include("utils.jl")

RNG = MersenneTwister(0)
PROC = ScatteringProcess((Electron(), Photon()), (Electron(), Photon(), Photon()))
MODEL = PerturbativeQED()
INPSL = FlatPhaseSpaceLayout(TwoBodyRestSystem())

include("impl/compton.jl")

@testset "Two Photon Compton" begin
    g = graph(PROC)

    # suppress type inference warnings; they don't matter here
    f = with_logger(ConsoleLogger(Logging.Error)) do
        get_compute_function(g, PROC, cpu_st(), @__MODULE__)
    end
    for i in 1:10
        input = gen_process_input(RNG, PROC)

        @test isapprox(f(input), two_compton_mat_el(input))
    end
end
