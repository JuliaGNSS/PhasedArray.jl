@testset "Pre Whitening filter" begin
    test_signal = complex.(randn(1000,2), randn(1000,2)) / sqrt(2)
    pre_whitening_filter = PhasedArray.calc_prewhitening_filter(test_signal)
    @test pre_whitening_filter ≈ I atol = 0.1
    @test filter(Matrix{Complex{Float64}}(I, 2, 2), test_signal) == test_signal
end

@testset "Amplitude filter" begin
    test_signal = complex.(randn(1000,2), randn(1000,2)) / sqrt(2)
    amplitude_filter = PhasedArray.calc_amplitude_filter(test_signal)
    @test amplitude_filter ≈ I atol = 0.1
end
