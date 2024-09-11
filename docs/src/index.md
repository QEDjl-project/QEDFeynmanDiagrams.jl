# QEDFeynmanDiagrams.jl

This project is dedicated to generating Feynman diagrams for scattering processes in perturbative QED. It is part of the [QEDjl-project](https://github.com/QEDjl-project) and relies on [QEDcore.jl](https://github.com/QEDjl-project/QEDcore.jl).
Furthermore, it can generate code to compute the matrix element for scattering processes, using [ComputableDAGs.jl](https://github.com/ComputableDAGs/ComputableDAGs.jl) and [RuntimeGeneratedFunctions.jl](https://github.com/SciML/RuntimeGeneratedFunctions.jl).

For example usage, see the [n-photon Compton](examples/compton.md) or [trident](examples/trident.md) pages.

The inner workings of the project are explained in the [manual](manual.md).
