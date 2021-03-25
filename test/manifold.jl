@testset "Manifold" begin
    f_0 = 1575420e3
    c₀ = 299_792_458
    λ = c₀ / f_0
    ant_pos = SVector(λ / 4 * SVector(1, 1, 0), λ / 4 * SVector(-1, 1, 0), λ / 4 * SVector(1, -1, 0), λ / 4 * SVector(-1, -1 , 0))
    manifold = @inferred IdealManifold(f_0, ant_pos)

    @test @inferred(get_num_ants(manifold)) == 4

    num_ants = length(EXAMPLE_LUT)
    manifold_real = RealManifold(EXAMPLE_LUT)
    @test @inferred(get_num_ants(manifold_real)) == num_ants
end