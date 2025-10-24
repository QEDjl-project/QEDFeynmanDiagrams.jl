using Random
using QEDbase.Mocks
using QEDcore
using ComputableDAGs
using QEDFeynmanDiagrams

using Logging

RNG = MersenneTwister(0)
PROC = bhabha
GRAPH = graph(PROC)
MODEL = MockModel()
INPSL = FlatPhaseSpaceLayout(TwoBodyRestSystem())
N = 128

# suppress type inference warnings; they don't matter here
f = with_logger(ConsoleLogger(Logging.Error)) do
    compute_function(GRAPH, PROC, cpu_st(), @__MODULE__)
end
k = eval(kernel(GRAPH, PROC, @__MODULE__))

@testset "Testing with $GPU_MODULE" for (GPU_MODULE, VECTOR_TYPE) in GPUS
    CDAG_GPU_T = GPU_TYPES_CDAG[GPU_MODULE]

    if isnothing(CDAG_GPU_T)
        @warn "$GPU_MODULE is not yet supported by ComputableDAGs.jl. Skipping GPU tests..."
        continue
    end

    @testset "Bhabha scattering on GPU ($MOM_EL_TYPE)" for MOM_EL_TYPE in
        GPU_FLOAT_TYPES[GPU_MODULE]
        input = [gen_process_input(RNG, PROC) for _ in 1:N]
        output = [zero(MOM_EL_TYPE) for _ in 1:N]
        gpu_input = VECTOR_TYPE(input)
        gpu_output = VECTOR_TYPE(output)

        expected_result = _ground_truth_bhabha.(input)

        @testset "generated kernel" begin
            k(get_backend(gpu_input), 32)(gpu_input, gpu_output; ndrange = N)

            @test eltype(gpu_output) == MOM_EL_TYPE
            @test isapprox(Vector(gpu_output), expected_result)
        end
    end
end
