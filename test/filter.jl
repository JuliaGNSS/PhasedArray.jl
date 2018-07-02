@testset "Pre Whitening filter" begin
    test_signal = complex.(randn(2,1000), randn(2,1000)) / sqrt(2)
    pre_whitening_filter = PhasedArray.calc_whitening_filter(test_signal)
    @test pre_whitening_filter ≈ eye(2) atol = 0.1
    @test filter(eye(2), test_signal) ==  test_signal
end