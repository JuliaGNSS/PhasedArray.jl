@testset "Pattern" begin
    manifold = IdealManifold(1575420e3, 0.1904 / 4 * SVector(SVector(1, 1, 0), SVector(-1, 1, 0), SVector(1, -1, 0), SVector(-1, -1, 0)))
    pattern = Pattern(manifold, norm, num_az = 360, num_el = 91, max_el = π / 2)
    @test pattern.azs == range(0, stop = 2π, length = 360)
    @test pattern.els == range(0, stop = π / 2, length = 91)
    @test pattern.values ≈ ones(91, 360) * 2
    @test pattern.max_el == π / 2
end

@testset "Pattern3D" begin
    manifold = IdealManifold(1575420e3, 0.1904 / 4 * SVector(SVector(1, 1, 0), SVector(-1, 1, 0), SVector(1, -1, 0), SVector(-1, -1, 0)))
    pattern = Pattern3D(manifold, norm, num_az = 360, num_el = 181, max_el = π)
    @test pattern.gains ≈ ones(181, 360) * 2
    @test pattern.max_el == 1π
end
