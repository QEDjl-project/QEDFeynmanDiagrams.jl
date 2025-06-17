using Test
using SafeTestsets

include("utils.jl")

# check if we run CPU tests (yes by default)
cpu_tests = _is_test_platform_active(["CI_QED_TEST_CPU", "TEST_CPU"], true)

if cpu_tests
    @safetestset "Number of diagrams" begin
        include("number_of_diagrams.jl")
    end

    @safetestset "Input Type" begin
        include("input_type.jl")
    end

    @safetestset "2-Photon Compton" begin
        include("two_photon_compton.jl")
    end

    @safetestset "Fermion Exchange" begin
        include("fermion_exchange.jl")
    end

    @safetestset "Synced Spins and Polarizations" begin
        include("synced_spin_pol.jl")
    end

    @safetestset "Test against Madgraph ground truths" begin
        include("madgraph/test_madgraph.jl")
    end
else
    @info "Skipping CPU tests"
end

begin
    @time @safetestset "GPU testing" begin
        include("gpu/runtests.jl")
    end
end
