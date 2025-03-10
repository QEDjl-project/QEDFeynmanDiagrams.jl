using SafeTestsets

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
