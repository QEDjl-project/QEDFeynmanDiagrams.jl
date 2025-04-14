function call_kernel(::Type{CUDAGPU}, k, in_vec, out_vec)
    t = 32
    b = length(in_vec) ÷ t

    CUDA.@sync (@cuda threads = t blocks = b k(in_vec, out_vec, length(in_vec)))
    return nothing
end
