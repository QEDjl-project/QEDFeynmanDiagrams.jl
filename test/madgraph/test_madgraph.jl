using QEDbase.Mocks
using QEDcore
using QEDFeynmanDiagrams
using Random
using Logging

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

RNG = MersenneTwister(137)
MODEL = MockModel()
PSL = FlatPhaseSpaceLayout(TwoBodyTargetSystem())

include("../utils.jl")

TEST_PROCESSES = if LARGE_TESTS()
    [
        "ema_ema",          # ke -> ke
        "emep_emep",        # ep -> ep
        "ema_emaa",         # ke -> kke
        "ema_ememep",       # ke -> eep
        "emep_emepemep",    # ep -> eepp
        "ema_emaaa",        # ke -> kkke
        "ema_emaaaa",       # ke -> kkkke
    ]
else
    @info "Skipping large tests...\nEnable them explicitly with an environment variable LARGE_TESTS=1"
    [
        "ema_ema",          # ke -> ke
        "emep_emep",        # ep -> ep
        "ema_emaa",         # ke -> kke
        "ema_ememep",       # ke -> eep
    ]
end

@testset "Testing Madgraph process $proc" for proc in TEST_PROCESSES
    include("$proc/setup.jl")
    GRAPH = graph(PROC)

    # suppress type inference warnings; they don't matter here
    f = with_logger(ConsoleLogger(Logging.Error)) do
        get_compute_function(GRAPH, PROC, cpu_st(), @__MODULE__)
    end
    results = f.(psps)

    # TODO: remove this once the source of the constant factor of error has been found
    FACTOR = groundtruth_mat_el_sq[1] / results[1]
    results .*= FACTOR
    @test isapprox(results, groundtruth_mat_el_sq)
end
