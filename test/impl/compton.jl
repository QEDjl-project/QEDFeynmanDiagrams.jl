using Combinatorics
using QEDFeynmanDiagrams

VERTEX = QEDFeynmanDiagrams.VERTEX

function assert_compton(proc::AbstractProcessDefinition)
    @assert number_particles(proc, Incoming(), Electron()) == 1 "there should be exactly one incoming electron"
    @assert number_particles(proc, Outgoing(), Electron()) == 1 "there should be exactly one outgoing electron"
    @assert number_particles(proc, Incoming(), Photon()) +
            number_particles(proc, Outgoing(), Photon()) +
            2 == number_particles(proc, Incoming()) + number_particles(proc, Outgoing()) "there should only be an incoming and outgoing electron, and photons"
end

function two_compton_diagram(
    in_el, out_el, ph_1, ph_2, ph_3, in_el_s, out_el_s, ph_1_p, ph_2_p, ph_3_p
)
    return base_state(Electron(), Outgoing(), -out_el, out_el_s) *
           (base_state(Photon(), Incoming(), ph_3, ph_3_p) * VERTEX) *
           propagator(Electron(), -out_el - ph_3) *
           (base_state(Photon(), Incoming(), ph_2, ph_2_p) * VERTEX) *
           propagator(Electron(), -out_el - ph_2 - ph_3) *
           (base_state(Photon(), Incoming(), ph_1, ph_1_p) * VERTEX) *
           base_state(Electron(), Incoming(), in_el, in_el_s)
end

function two_compton_spin_pol_comb(psp, spin_pols)
    in_el = momentum(psp, Incoming(), Electron())
    in_el_s = spin_pols[1][1]
    out_el = -momentum(psp, Outgoing(), Electron())
    out_el_s = spin_pols[2][1]

    ph_1 = momentum(psp, Incoming(), Photon(), Val(1))
    ph_2 = -momentum(psp, Outgoing(), Photon(), Val(1))
    ph_3 = -momentum(psp, Outgoing(), Photon(), Val(2))

    ph_1_p = spin_pols[1][2]
    ph_3_p = spin_pols[2][2]
    ph_2_p = spin_pols[2][3]

    diag1 = two_compton_diagram(
        in_el, out_el, ph_1, ph_2, ph_3, in_el_s, out_el_s, ph_1_p, ph_2_p, ph_3_p
    )
    diag2 = two_compton_diagram(
        in_el, out_el, ph_1, ph_3, ph_2, in_el_s, out_el_s, ph_1_p, ph_3_p, ph_2_p
    )
    diag3 = two_compton_diagram(
        in_el, out_el, ph_2, ph_1, ph_3, in_el_s, out_el_s, ph_2_p, ph_1_p, ph_3_p
    )
    diag4 = two_compton_diagram(
        in_el, out_el, ph_2, ph_3, ph_1, in_el_s, out_el_s, ph_2_p, ph_3_p, ph_1_p
    )
    diag5 = two_compton_diagram(
        in_el, out_el, ph_3, ph_1, ph_2, in_el_s, out_el_s, ph_3_p, ph_1_p, ph_2_p
    )
    diag6 = two_compton_diagram(
        in_el, out_el, ph_3, ph_2, ph_1, in_el_s, out_el_s, ph_3_p, ph_2_p, ph_1_p
    )

    return abs2(diag1 + diag2 + diag3 + diag4 + diag5 + diag6)
end

function two_compton_mat_el(psp::AbstractPhaseSpacePoint)
    proc = process(psp)
    assert_compton(proc)
    @assert number_particles(proc, Outgoing(), Photon()) == 2 "there should be exactly two outgoing photons"
    @assert number_particles(proc, Incoming(), Photon()) == 1 "there should be exactly one incoming photon"
    @assert number_incoming_particles(proc) == 2 "only 2 incoming particles are allowed"
    @assert incoming_particles(proc)[1] == Electron() "the first incoming particle should be the electron"
    @assert outgoing_particles(proc)[1] == Electron() "the first outgoing particle should be the electron"

    TYPE = Float64 #momentum_eltype(psp)

    cum_sum = zero(TYPE)
    for sp_comb in spin_pols_iter(proc)
        @debug sp_comb

        cum_sum += two_compton_spin_pol_comb(psp, sp_comb)
    end

    return cum_sum
end
