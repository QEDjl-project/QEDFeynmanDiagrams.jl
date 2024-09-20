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

# make sure we're not keeping all these graphs in memory
GC.gc()

@testset "Trident-like process with $n produced pairs" for n in (1, 2)
    # sync the positrons and the electrons with each other, across input and output
    proc_synced = ScatteringProcess(
        (Electron(), Photon()),
        (Electron(), ntuple(_ -> Electron(), n)..., ntuple(_ -> Positron(), n)...),
        (SyncedSpin(1), PolX()),
        (SyncedSpin(1), ntuple(_ -> SyncedSpin(1), n)..., ntuple(_ -> SyncedSpin(2), n)...),
    )

    # this gives 4 combinations: up-up, up-down, down-up, and down-down
    procs = ScatteringProcess[]
    for (spin1, spin2) in Iterators.product((SpinUp(), SpinDown()), (SpinUp(), SpinDown()))
        push!(
            procs,
            ScatteringProcess(
                (Electron(), Photon()),
                (Electron(), ntuple(_ -> Electron(), n)..., ntuple(_ -> Positron(), n)...),
                (spin1, PolX()),
                (spin1, ntuple(_ -> spin1, n)..., ntuple(_ -> spin2, n)...),
            ),
        )
    end

    g_synced = generate_DAG(proc_synced)
    graphs = generate_DAG.(procs)

    for g in graphs[2:end]
        @test length(graphs[1].nodes) == length(g.nodes)
    end

    f_synced = get_compute_function(g_synced, proc_synced, cpu_st(), @__MODULE__)
    functions = get_compute_function.(graphs, procs, Ref(cpu_st()), Ref(@__MODULE__))

    inputs_synced = [gen_process_input(RNG, proc_synced) for _ in 1:1000]

    # have to cast these for their respective processes
    inputs_unsynced = Vector{PhaseSpacePoint}[]
    for p in procs
        push!(
            inputs_unsynced,
            [
                PhaseSpacePoint(
                    p,
                    model(psp),
                    phase_space_definition(psp),
                    momenta(psp, Incoming()),
                    momenta(psp, Outgoing()),
                ) for psp in inputs_synced
            ],
        )
    end

    results_synced = f_synced.(inputs_synced)
    results_unsynced = Vector{Float64}[]
    for (f, inputs) in zip(functions, inputs_unsynced)
        push!(results_unsynced, f.(inputs))
    end

    results_summed = [sum(args) for args in zip(results_unsynced...)]

    @test isapprox(results_synced, results_summed)
end

GC.gc()
