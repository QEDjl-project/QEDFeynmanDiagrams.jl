# # n-Photon Compton Scattering Process

# In this file, we set up an n-photon Compton scattering process. A Compton scattering
# process looks like $k^n e^- \to k e^-$.

# You can download this file as a [jupyter notebook](compton.ipynb).

using QEDFeynmanDiagrams

# We need some of the packages of the [QEDjl-project](https://github.com/QEDjl-project) for base
# functionality and the `ScatteringProcess` type.
using QEDcore
using QEDprocesses

# Let's decide how many photons our electron interacts with:
n = 4;

# Now we setup the scattering process accordingly. We consider all spin/polarization
# combinations of the particles except for the incoming photons, where the polarizations are synced using [`QEDbase.SyncedPolarization`](@extref). 
# This emulates all synced photons having the same, but still indefinite, polarization, for example from a laser.
proc = ScatteringProcess(
    (Electron(), ntuple(_ -> Photon(), n)...),     # incoming particles
    (Electron(), Photon()),                        # outgoing particles
    (AllSpin(), ntuple(_ -> SyncedPol(1), n)...),  # incoming particle spin/pols
    (AllSpin(), AllPol()),                         # outgoing particle spin/pols
)

# The [`feynman_diagrams`](@ref) function returns an iterator for all possible Feynman diagrams
# for this scattering process. With its `length` overload, we can check how many diagrams
# there are. For an n-photon Compton process with `n` incoming photons, this should be
# $(n+1)!$.
length(feynman_diagrams(proc))

# Next, we can generate the DAG representing the computation for our scattering process'
# squared matrix element. This uses [`ComputableDAGs.jl`](https://github.com/ComputableDAGs/ComputableDAGs.jl).
dag = generate_DAG(proc)

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

# With the DAG, the process, and `RuntimeGeneratedFunctions` initalized,
# we can now generate the actual computable function:
func = get_compute_function(dag, proc, cpu_st(), @__MODULE__);

# Now we need an input for the function, which is a [`QEDcore.PhaseSpacePoint`](@extref).
# For now, we generate random momenta for every particle. In the future, QEDevents
# will be able to generate physical `PhaseSpacePoint`s.
psp = PhaseSpacePoint(
    proc,
    PerturbativeQED(),
    PhasespaceDefinition(SphericalCoordinateSystem(), ElectronRestFrame()),
    tuple((rand(SFourMomentum) for _ in 1:number_incoming_particles(proc))...),
    tuple((rand(SFourMomentum) for _ in 1:number_outgoing_particles(proc))...),
)

# Finally, we can test that the function actually runs and computes something by
# simply calling it on the `PhaseSpacePoint`:
func(psp)

# If we want, we can benchmark the execution speed too:
using BenchmarkTools
@benchmark func($psp)
