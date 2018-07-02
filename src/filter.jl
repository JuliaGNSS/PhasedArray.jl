function calc_whitening_filter(signal)
    Rxx = signal * ctranspose(signal) / size(signal, 2)
    return Rxx^-0.5
end

function filter(filter_matrix::AbstractArray, signal::Array{Complex{Float64},2})
    filter_matrix * signal
end