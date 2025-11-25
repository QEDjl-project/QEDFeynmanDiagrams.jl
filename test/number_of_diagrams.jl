using QEDbase.Mocks
using QEDFeynmanDiagrams

include("utils.jl")

RNG = MersenneTwister(0)

@testset "Compton-like process with $n incoming photons" for n in 1:7
    proc = MockProcess((Electron(), ntuple(_ -> Photon(), n)...), (Electron(), Photon()))

    @test factorial(n + 1) == number_of_diagrams(proc)
end

@testset "Trident-like process with $n produced pairs" for n in 1:7
    proc = MockProcess(
        (Electron(), ntuple(_ -> Photon(), n)...),
        (Electron(), Electron(), Positron(), Photon()),
    )

    @test factorial(n + 4, 3) * 2 == number_of_diagrams(proc)
end

@testset "Invalid Processes" begin
    @test number_of_diagrams(MockProcess((Electron(), Electron()), (Positron(), Electron()))) == 0
    @test number_of_diagrams(MockProcess((Photon(), Photon()), (Photon(), Photon()))) == 0 # diagrams exist, but only with loops
    @test number_of_diagrams(MockProcess((Photon(), Photon()), (Photon(), Positron()))) == 0
end
