@testset "Pattern" begin
    manifold = IdealManifold(1575420e3, 0.1904 / 4 * SVector(SVector(1, 1, 0), SVector(-1, 1, 0), SVector(1, -1, 0), SVector(-1, -1, 0)))
    pattern = Pattern(manifold, norm, num_az = 360, num_el = 91, max_el = π / 2)
    @test first(pattern.azs) == 0.0
    @test step(pattern.azs) == 2π/360
    @test length(pattern.azs) == 360
    @test isa(pattern.azs, AbstractRange)
    @test first(pattern.els) == 0.0
    @test step(pattern.els) == π/2 / 90
    @test length(pattern.els) == 91
    @test isa(pattern.els, AbstractRange)
    @test pattern.values ≈ ones(91, 360) * 2
    @test pattern.max_el == π / 2
end

@testset "Pattern3D" begin
    manifold = IdealManifold(1575420e3, 0.1904 / 4 * SVector(SVector(1, 1, 0), SVector(-1, 1, 0), SVector(1, -1, 0), SVector(-1, -1, 0)))
    pattern = Pattern3D(manifold, norm, num_az = 360, num_el = 181, max_el = π)
    @test pattern.gains ≈ ones(181, 360) * 2
    @test pattern.max_el == 1π
end
