function calc_prewhitening_filter(Rxx)
    F = eigen(Rxx)
    F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors' .* sqrt(mean_power(F.values))
end

function calc_amplitude_filter(Rxx)
    F = eigen(Rxx)
    F.vectors * Diagonal(1 ./ F.values) * F.vectors' .* mean_power(F.values)
end

function calc_eigen_beamformer(Rxx)
    eigvecs(Rxx)[:,end]
end

function calc_variance_covariance(signal)
    Hermitian(signal' * signal) / size(signal, 1)
end

function mean_power(eigenvalues)
    1 / (sum(1 ./ eigenvalues) / length(eigenvalues))
end
