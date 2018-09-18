function calc_prewhitening_filter(signal)
    Rxx = signal' * signal / size(signal, 1)
    Rxx^-0.5
end

function calc_amplitude_filter(signal)
    Rxx = signal' * signal / size(signal, 1)
    Rxx^-1
end

function filter(filter_matrix::AbstractArray, signal::Array{Complex{Float64},2})
    signal * filter_matrix
end
