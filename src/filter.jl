function calc_prewhitening_filter(signal)
    Rxx = calc_variance_covariance(signal)
    F = eigen(Rxx)
    F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors' .* sqrt(mean_power(F.values))
end

function calc_amplitude_filter(signal)
    Rxx = calc_variance_covariance(signal)
    F = eigen(Rxx)
    F.vectors * Diagonal(1 ./ F.values) * F.vectors' .* mean_power(F.values)
end

function filter(filter_matrix::AbstractArray, signal::Array{Complex{Float64},2})
    signal * filter_matrix'
end

function calc_eigen_beamformer(signal)
    Rxx = calc_variance_covariance(signal)
    eigvecs(Rxx)[:,end]
end

function calc_variance_covariance(signal)
    Hermitian(signal' * signal) / size(signal, 1)
end

function mean_power(eigenvalues)
    1 / (sum(1 ./ eigenvalues) / length(eigenvalues))
end
