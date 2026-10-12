function calc_covariance(A, signal_powers, noise_power)
    Hermitian(A * Diagonal(signal_powers) * A' + I * noise_power)
end

function samv2(
    manifold::AbstractManifold,
    covariance_matrix::Hermitian;
    num_doas,
    min_phi = -10 * π / 180,
    max_iterations = 100,
    break_threshold = 1e-3
)
    doas = create_points_on_sphere(num_doas, min_phi)
    steer_vecs = map(doa -> get_steer_vec(manifold, doa), doas)
    A = reduce(hcat, steer_vecs)
    signal_powers = map(steer_vec -> abs(steer_vec' * covariance_matrix * steer_vec) / norm(steer_vec)^4, steer_vecs)
    noise_power = first(eigvals(covariance_matrix))
    for i = 1:max_iterations
        estimated_covariance_matrix = calc_covariance(A, signal_powers, noise_power)
        next_signal_powers = map(signal_powers, steer_vecs) do signal_power, steer_vec
            abs(signal_power * steer_vec' * (estimated_covariance_matrix \ (covariance_matrix * (estimated_covariance_matrix \ steer_vec))) /
                (steer_vec' * (estimated_covariance_matrix \ steer_vec)))
        end
        noise_power = abs(tr(estimated_covariance_matrix \ (estimated_covariance_matrix \ covariance_matrix)) /
            tr(inv(estimated_covariance_matrix)^2))
        normed_diff = norm(signal_powers - next_signal_powers) / norm(signal_powers)
        signal_powers = next_signal_powers
        normed_diff < break_threshold && break;
    end
    doas, signal_powers
end