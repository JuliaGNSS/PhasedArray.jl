function calc_whitening_filter(signal)
    Rxx = signal' * signal / size(signal, 1)
    Rxx^-0.5
end

function filter(filter_matrix::AbstractArray, signal::Array{Complex{Float64},2})
    signal * filter_matrix
end
