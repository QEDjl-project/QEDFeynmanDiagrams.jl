# file for testing that the generated input_type of a generated dag is correct

using Random
using QEDcore
using QEDprocesses
using ComputableDAGs
using QEDFeynmanDiagrams

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

include("utils.jl")

RNG = MersenneTwister(0)

@testset "Compton-like process with $n incoming photons" for n in (1, 2, 3, 4)
    proc = ScatteringProcess(
        (Electron(), ntuple(_ -> Photon(), n)...),
        (Electron(), Photon()),
        (AllSpin(), ntuple(_ -> PolX(), n)...),
        (AllSpin(), AllPol()),
    )

    for n_other in (1, 2, 3, 4)
        n_other_proc = ScatteringProcess(
            (Electron(), ntuple(_ -> Photon(), n_other)...),
            (Electron(), Photon()),
            (AllSpin(), ntuple(_ -> PolX(), n_other)...),
            (AllSpin(), AllPol()),
        )
        input = gen_process_input(RNG, n_other_proc)

        if n_other == n
            @test typeof(input) <: input_type(proc)
        else
            @test !(typeof(input) <: input_type(proc))
        end
    end
end
