# file for testing that synced spins and polarizations are handled correctly

using Random
using QEDcore
using QEDprocesses
using ComputableDAGs
using QEDFeynmanDiagrams

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

include("utils.jl")

RNG = MersenneTwister(0)

@testset "Compton-like process with $n incoming photons" for n in (1, 2, 3)
    proc_synced = ScatteringProcess(
        (Electron(), ntuple(_ -> Photon(), n)...),
        (Electron(), Photon()),
        (AllSpin(), ntuple(_ -> SyncedPol(1), n)...),
        (AllSpin(), AllPol()),
    )

    proc_polx = ScatteringProcess(
        (Electron(), ntuple(_ -> Photon(), n)...),
        (Electron(), Photon()),
        (AllSpin(), ntuple(_ -> PolX(), n)...),
        (AllSpin(), AllPol()),
    )
    proc_poly = ScatteringProcess(
        (Electron(), ntuple(_ -> Photon(), n)...),
        (Electron(), Photon()),
        (AllSpin(), ntuple(_ -> PolY(), n)...),
        (AllSpin(), AllPol()),
    )

    g_synced = generate_DAG(proc_synced)
    g_polx = generate_DAG(proc_polx)
    g_poly = generate_DAG(proc_poly)

    @test length(g_polx.nodes) == length(g_poly.nodes)
    @test length(g_polx.nodes) < length(g_synced.nodes) < 2 * length(g_polx.nodes)

    f_synced = get_compute_function(g_synced, proc_synced, cpu_st(), @__MODULE__)
    f_polx = get_compute_function(g_polx, proc_polx, cpu_st(), @__MODULE__)
    f_poly = get_compute_function(g_poly, proc_poly, cpu_st(), @__MODULE__)

    inputs_synced = [gen_process_input(RNG, proc_synced) for _ in 1:1000]
    # have to cast these for their respective processes
    inputs_polx = [
        PhaseSpacePoint(
            proc_polx,
            model(psp),
            phase_space_definition(psp),
            momenta(psp, Incoming()),
            momenta(psp, Outgoing()),
        ) for psp in inputs_synced
    ]
    inputs_poly = [
        PhaseSpacePoint(
            proc_poly,
            model(psp),
            phase_space_definition(psp),
            momenta(psp, Incoming()),
            momenta(psp, Outgoing()),
        ) for psp in inputs_synced
    ]

    results_synced = f_synced.(inputs_synced)
    results_polx = f_polx.(inputs_polx)
    results_poly = f_poly.(inputs_poly)

    @test isapprox(results_synced, results_polx .+ results_poly)
end
