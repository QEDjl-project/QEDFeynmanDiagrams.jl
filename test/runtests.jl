using SafeTestsets

@safetestset "Input Type" begin
    include("input_type.jl")
end

@safetestset "Fermion Exchange" begin
    include("fermion_exchange.jl")
end

@safetestset "Synced Spins and Polarizations" begin
    include("synced_spin_pol.jl")
end
