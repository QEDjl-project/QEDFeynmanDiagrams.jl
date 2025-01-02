using QEDprocesses
using QEDFeynmanDiagrams

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

include("utils.jl")

RNG = MersenneTwister(0)

@testset "Compton-like process with $n incoming photons" for n in 1:7
    proc = ScatteringProcess(
        (Electron(), ntuple(_ -> Photon(), n)...), (Electron(), Photon())
    )

    @test factorial(n + 1) == number_of_diagrams(proc)
end

@testset "Trident-like process with $n produced pairs" for n in 1:7
    proc = ScatteringProcess(
        (Electron(), ntuple(_ -> Photon(), n)...),
        (Electron(), Electron(), Positron(), Photon()),
    )

    @test factorial(n + 4, 3) * 2 == number_of_diagrams(proc)
end
