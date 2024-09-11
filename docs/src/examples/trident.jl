# # n-Pair Trident Scattering Process

# In this file, we set up an n-pair trident scattering process. A trident
# process looks like $k e^- \to e^- (e^- e^+)^n$.

# You can download this file as a [jupyter notebook](trident.ipynb).

using QEDFeynmanDiagrams

# We need some of the packages of the [QEDjl-project](https://github.com/QEDjl-project) for base
# functionality and the `ScatteringProcess` type.
using QEDcore
using QEDprocesses

# Let's decide how many pairs our trident should produce:
n = 2;

# Now we setup the scattering process accordingly. We only consider a single spin/polarization
# combination here. For an example with many spin and polarization combinations, refer to the
# [n-photon Compton example](compton.md)
proc = ScatteringProcess(
    (Electron(), Photon()),                                                         # incoming particles
    (Electron(), ntuple(_ -> Electron(), n)..., ntuple(_ -> Positron(), n)...),     # outgoing particles
    (SpinUp(), PolX()),                                                             # incoming particle spin/pols
    (SpinUp(), ntuple(_ -> SpinUp(), 2 * n)...),                                    # outgoing particle spin/pols
)

# The [`feynman_diagrams`](@ref) function returns an iterator for all possible Feynman diagrams
# for this scattering process. With its `length` overload, we can check how many diagrams
# there are.
length(feynman_diagrams(proc))

# Next, we can generate the DAG representing the computation for our scattering process'
# squared matrix element. This uses `ComputableDAGs.jl`.
dag = generate_DAG(proc)

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
