@testset "Pre Whitening filter" begin
    test_signal = complex.(randn(1000, 2), randn(1000, 2)) / sqrt(2)
    Rxx = calc_variance_covariance(test_signal)
    pre_whitening_filter = calc_prewhitening_filter(Rxx)
    @test pre_whitening_filter ≈ I atol = 0.1

    test_signal = complex.(randn(1000, 2), randn(1000, 2)) / sqrt(2) * sqrt(100) # Power: 100
    Rxx = calc_variance_covariance(test_signal)
    pre_whitening_filter = calc_prewhitening_filter(Rxx)
    @test pre_whitening_filter ≈ I atol = 0.1
end

@testset "Amplitude filter" begin
    test_signal = complex.(randn(1000, 2), randn(1000, 2)) / sqrt(2)
    Rxx = calc_variance_covariance(test_signal)
    amplitude_filter = calc_amplitude_filter(Rxx)
    @test amplitude_filter ≈ I atol = 0.1

    test_signal = complex.(randn(1000, 2), randn(1000, 2)) / sqrt(2) * sqrt(100) # Power: 100
    Rxx = calc_variance_covariance(test_signal)
    pre_whitening_filter = calc_amplitude_filter(Rxx)
    @test pre_whitening_filter ≈ I atol = 0.1
end

@testset "Eigen beamformer" begin
    test_signal = complex.(randn(1000, 4), randn(1000, 4)) / sqrt(2) / 5 .+ ones(1000, 4)
    Rxx = calc_variance_covariance(test_signal)
    beamformer = calc_eigen_beamformer(Rxx)
    @test beamformer ≈ -1 * [0.5, 0.5, 0.5, 0.5] atol = 0.03
end
