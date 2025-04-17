function call_kernel(::Type{ROCmGPU}, k, in_vec, out_vec)
    t = 32
    b = length(in_vec) ÷ t

    AMDGPU.@sync (@roc groupsize = t gridsize = b k(in_vec, out_vec, length(in_vec)))
    return nothing
end
