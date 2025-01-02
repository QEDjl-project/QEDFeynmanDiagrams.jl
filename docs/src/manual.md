# Manual

## Usage

!!! note 
    This project uses [ComputableDAGs.jl](https://github.com/ComputableDAGs/ComputableDAGs.jl), which uses [RuntimeGeneratedFunctions.jl](https://github.com/SciML/RuntimeGeneratedFunctions.jl). The latter requires an initialization step, so when using this project, make sure to have
    ```julia
    using RuntimeGeneratedFunctions
    RuntimeGeneratedFunctions.init(@__MODULE__)
    ```
    somewhere in the preamble of your code or at the beginning of your REPL session. Otherwise the function generation will not work.

This package currently exports two functions, [`graph`](@ref) and [`number_of_diagrams`](@ref). Both of these take a [`QEDbase.AbstractProcessDefinition`](@extref) as input. `graph` returns a [`ComputableDAGs.DAG`](@extref) representing the computation of the squared matrix element for [`QEDcore.PhaseSpacePoint`](@extref)s. `number_of_diagrams` is a utility function that returns the number of valid tree-level Feynman diagrams for the given process.

For usage examples, please refer to the examples section of the docs.
