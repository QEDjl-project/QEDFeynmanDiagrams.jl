# # n-Photon Compton Scattering Process

# In this file, we set up an n-photon Compton scattering process. A Compton scattering
# process looks like $k^n e^- \to k e^-$.

# You can download this file as a [jupyter notebook](compton.ipynb).

using QEDFeynmanDiagrams

# We need QEDcore of the [QEDjl-project](https://github.com/QEDjl-project) for base
# functionality and a process type, for which we can use the `Mocks` submodule
# from QEDbase for this tutorial. Downstream, a `ScatteringProcess` from QEDprocesses.jl
# could be used, for example.
using QEDcore
using QEDbase.Mocks

# Let's decide how many photons our electron interacts with:
n = 4;

# Now we setup the scattering process accordingly. We consider all spin/polarization
# combinations of the particles except for the incoming photons, where the polarizations are synced using [`QEDbase.SyncedPolarization`](@extref). 
# This emulates all synced photons having the same, but still indefinite, polarization, for example from a laser.
# !!! note
#     Currently, this process uses outgoing photons instead of incoming photons, because there is not yet a
#     `PhaseSpaceLayout` for more than two incoming particles in QEDcore.jl. See issue https://github.com/QEDjl-project/QEDcore.jl/issues/103
proc = QEDbase.Mocks.MockProcessSP(
    (Electron(), Photon()),                        # incoming particles
    (Electron(), ntuple(_ -> Photon(), n)...),     # outgoing particles
    (AllSpin(), AllPol()),                         # incoming particle spin/pols
    (AllSpin(), ntuple(_ -> SyncedPol(1), n)...),  # outgoing particle spin/pols
)

# The [`number_of_diagrams`](@ref) function returns how many diagrams there are for a given process.
# For an n-photon Compton process with `n` incoming photons, this should be $(n+1)!$.
number_of_diagrams(proc)

# Next, we can generate the DAG representing the computation for our scattering process'
# squared matrix element. This uses [`ComputableDAGs.jl`](https://github.com/ComputableDAGs/ComputableDAGs.jl).
dag = graph(proc)

# In this graph output you can see the number of nodes necessary to compute.
# Note that for larger processes, the number of total nodes can be *lower* than
# the number of Feynman diagrams, even with the added complexity of considering
# multiple spin and polarization combinations. This is the result of efficient
# reuse of reappearing parts of Feynman diagrams.

# To continue, we will need [`ComputableDAGs.jl`](https://github.com/ComputableDAGs/ComputableDAGs.jl). Since `ComputableDAGs.jl` uses 
# `RuntimeGeneratedFunction`s as the return type of [`ComputableDAGs.get_compute_function`](@extref), we need
# to initialize it in our current module.
using ComputableDAGs
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# Now we need an input for the function, which is a [`QEDcore.PhaseSpacePoint`](@extref).
# For now, we generate random momenta for every particle. In the future, QEDevents
# will be able to generate physical `PhaseSpacePoint`s.
psp = PhaseSpacePoint(
    proc,
    MockModel(),
    FlatPhaseSpaceLayout(TwoBodyRestSystem()),
    tuple((rand(SFourMomentum) for _ in 1:number_incoming_particles(proc))...),
    tuple((rand(SFourMomentum) for _ in 1:number_outgoing_particles(proc))...),
)

# With the DAG, the process, `RuntimeGeneratedFunctions` initialized, and an input type to use,
# we can now generate the actual computable function:
func = get_compute_function(
    dag, proc, cpu_st(), @__MODULE__; concrete_input_type=typeof(psp)
);

# Finally, we can test that the function actually runs and computes something by
# simply calling it on the `PhaseSpacePoint`:
func(psp)

# We can benchmark the execution speed too:
using BenchmarkTools
@benchmark func($psp)
