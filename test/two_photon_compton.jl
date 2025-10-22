# file for testing that fermion exchange negation is handled correctly

using Random
using QEDbase.Mocks
using QEDcore
using ComputableDAGs
using QEDFeynmanDiagrams

using Logging

include("utils.jl")

RNG = MersenneTwister(0)
PROC = MockProcess((Electron(), Photon()), (Electron(), Photon(), Photon()))
MODEL = MockModel()
INPSL = FlatPhaseSpaceLayout(TwoBodyRestSystem())

include("impl/compton.jl")

@testset "Two Photon Compton" begin
    g = graph(PROC)

    # suppress type inference warnings; they don't matter here
    f = with_logger(ConsoleLogger(Logging.Error)) do
        compute_function(g, PROC, cpu_st(), @__MODULE__)
    end
    for i in 1:10
        input = gen_process_input(RNG, PROC)

        @test isapprox(f(input), two_compton_mat_el(input))
    end
end
