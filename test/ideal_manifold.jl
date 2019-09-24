@testset "Ideal manifold" begin
    f_0 = 1575420e3
    c₀ = 299_792_458
    λ = c₀ / f_0
    ant_pos = SVector(λ / 4 * SVector(1, 1, 0), λ / 4 * SVector(-1, 1, 0), λ / 4 * SVector(1, -1, 0), λ / 4 * SVector(-1, -1 , 0))
    manifold = @inferred IdealManifold(ant_pos, f_0)

    @test @inferred(get_steer_vec(manifold, SVector(0,0,1))) ≈ [1, 1, 1, 1]
    @test @inferred(get_steer_vec(manifold, SVector(0,1,0))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0))) ≈ [1im, -1im, 1im, -1im]

    @test @inferred(get_steer_vec(manifold, Spherical(1.0, 0.0, 0.0))) ≈ [1im, -1im, 1im, -1im]
    @test @inferred(get_steer_vec(manifold, Spherical(1.0, 0.0, 0.0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]

    manifold = @inferred IdealManifold(SVector{1, SVector{3, Float64}}(SVector(0.0,0.0,0.0)), f_0)
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0))) == 1.0
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0), RotXYZ(0.0,0.0,0.0))) == 1.0
end
