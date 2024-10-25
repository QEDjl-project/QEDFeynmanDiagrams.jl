# file for testing that fermion exchange negation is handled correctly

using Random
using QEDcore
using QEDprocesses
using ComputableDAGs
using QEDFeynmanDiagrams

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

VERTEX = QEDFeynmanDiagrams.VERTEX

include("utils.jl")

RNG = MersenneTwister(0)

bhabha = ScatteringProcess(
    (Electron(), Positron()),
    (Electron(), Positron()),
    (AllSpin(), AllSpin()),
    (AllSpin(), AllSpin()),
)

BHABHA = typeof(bhabha)

# groundtruth for ep -> ep
function _ground_truth_bhabha(psp::PhaseSpacePoint{BHABHA,PerturbativeQED})
    e_i = psp[Incoming(), 1]
    e_o = psp[Outgoing(), 1]
    p_i = psp[Incoming(), 2]
    p_o = psp[Outgoing(), 2]

    # propagators
    prop_ei_pi = propagator(Photon(), momentum(e_i) + momentum(p_i))
    prop_ei_eo = propagator(Photon(), momentum(e_i) - momentum(e_o))

    sum = 0.0

    for ((e_i_s, p_i_s), (e_o_s, p_o_s)) in spin_pols_iter(process(psp))
        # base states
        e_i_b = base_state(
            particle_species(e_i), particle_direction(e_i), momentum(e_i), e_i_s
        )
        e_o_b = base_state(
            particle_species(e_o), particle_direction(e_o), momentum(e_o), e_o_s
        )
        p_i_b = base_state(
            particle_species(p_i), particle_direction(p_i), momentum(p_i), p_i_s
        )
        p_o_b = base_state(
            particle_species(p_o), particle_direction(p_o), momentum(p_o), p_o_s
        )

        # diagram 1: e_i with p_i - photon - e_o with p_o
        diagram1 = (p_i_b * VERTEX * e_i_b) * prop_ei_pi * (e_o_b * VERTEX * p_o_b)
        diagram2 = (e_o_b * VERTEX * e_i_b) * prop_ei_eo * (p_i_b * VERTEX * p_o_b)
        sum += abs2(diagram1 - diagram2)
    end

    return sum
end

@testset "Bhabha Scattering ep -> ep" begin
    g = graph(bhabha)
    f = get_compute_function(g, bhabha, cpu_st(), @__MODULE__)
    for i in 1:10
        input = gen_process_input(RNG, bhabha)

        # println("gt1 : $(_ground_truth_bhabha(input))")
        # println("gt2 : $(_ground_truth_bhabha_test(input))\n")

        @test isapprox(f(input), _ground_truth_bhabha(input))
    end
end
