using Combinatorics
using QEDFeynmanDiagrams
using QEDcore

using QEDFeynmanDiagrams: VERTEX

bhabha = Mocks.MockProcessSP(
    (Electron(), Positron()),
    (Electron(), Positron()),
    (AllSpin(), AllSpin()),
    (AllSpin(), AllSpin()),
)

BHABHA = typeof(bhabha)

function assert_bhabha(proc::AbstractProcessDefinition)
    @assert number_particles(proc, Incoming(), Electron()) == 1 "there should be exactly one incoming electron"
    @assert number_particles(proc, Outgoing(), Electron()) == 1 "there should be exactly one outgoing electron"
    @assert number_particles(proc, Incoming(), Positron()) == 1 "there should be exactly one incoming positron"
    @assert number_particles(proc, Outgoing(), Positron()) == 1 "there should be exactly one outgoing positron"
    @assert number_particles(proc, Incoming(), Photon()) == 0 "there should be zero incoming photons"
    @assert number_particles(proc, Outgoing(), Photon()) == 0 "there should be zero outgoing photons"
    return nothing
end

# ground truth for ep -> ep
function _ground_truth_bhabha(psp::PhaseSpacePoint)
    assert_bhabha(process(psp))

    e_i = psp[Incoming(), 1]
    e_o = psp[Outgoing(), 1]
    p_i = psp[Incoming(), 2]
    p_o = psp[Outgoing(), 2]

    T = momentum_eltype(psp)

    # propagators
    prop_ei_pi = propagator(Photon(), momentum(e_i) + momentum(p_i))
    prop_ei_eo = propagator(Photon(), momentum(e_i) - momentum(e_o))

    sum = zero(T)

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
        diagram1 = (p_i_b * VERTEX(T) * e_i_b) * prop_ei_pi * (e_o_b * VERTEX(T) * p_o_b)
        diagram2 = (e_o_b * VERTEX(T) * e_i_b) * prop_ei_eo * (p_i_b * VERTEX(T) * p_o_b)
        sum += abs2(diagram1 - diagram2)
    end

    return sum
end
