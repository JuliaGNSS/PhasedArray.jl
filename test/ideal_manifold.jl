@testset "Ideal manifold" begin
    f_0 = 1575420e3
    c₀ = 299_792_458
    λ = c₀ / f_0
    ant_pos = SVector(λ / 4 * SVector(1, 1, 0), λ / 4 * SVector(-1, 1, 0), λ / 4 * SVector(1, -1, 0), λ / 4 * SVector(-1, -1 , 0))
    manifold = @inferred IdealManifold(f_0, ant_pos)

    @test typeof(manifold) <: AbstractManifold{4}
    @test @inferred(get_steer_vec(manifold, SVector(0,0,1))) ≈ [1, 1, 1, 1]
    @test @inferred(get_steer_vec(manifold, SVector(0,1,0))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0))) ≈ [1im, -1im, 1im, -1im]

    @test @inferred(get_steer_vec(manifold, Spherical(1.0, 0.0, 0.0))) ≈ [1im, -1im, 1im, -1im]
    @test @inferred(get_steer_vec(manifold, Spherical(1.0, 0.0, 0.0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]

    manifold = @inferred IdealManifold(f_0, λ / 4 * SVector(1, 1, 0), λ / 4 * SVector(-1, 1, 0), λ / 4 * SVector(1, -1, 0), λ / 4 * SVector(-1, -1 , 0))

    @test typeof(manifold) <: AbstractManifold{4}
    @test @inferred(get_steer_vec(manifold, SVector(0,0,1))) ≈ [1, 1, 1, 1]
    @test @inferred(get_steer_vec(manifold, SVector(0,1,0))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0))) ≈ [1im, -1im, 1im, -1im]

    @test @inferred(get_steer_vec(manifold, Spherical(1.0, 0.0, 0.0))) ≈ [1im, -1im, 1im, -1im]
    @test @inferred(get_steer_vec(manifold, Spherical(1.0, 0.0, 0.0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]

    manifold = @inferred IdealManifold(f_0, λ / 4 * @SMatrix([1 -1 1 -1; 1 1 -1 -1; 0 0 0 0]))

    @test typeof(manifold) <: AbstractManifold{4}
    @test @inferred(get_steer_vec(manifold, SVector(0,0,1))) ≈ [1, 1, 1, 1]
    @test @inferred(get_steer_vec(manifold, SVector(0,1,0))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0))) ≈ [1im, -1im, 1im, -1im]

    @test @inferred(get_steer_vec(manifold, Spherical(1.0, 0.0, 0.0))) ≈ [1im, -1im, 1im, -1im]
    @test @inferred(get_steer_vec(manifold, Spherical(1.0, 0.0, 0.0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]

    manifold = @inferred IdealManifold()
    @test typeof(manifold) <: AbstractManifold{1}
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0))) == 1.0
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0), RotXYZ(0.0,0.0,0.0))) == 1.0

    num_ant_x = 2
    dist_x = 0.5
    manifold = IdealManifold(1, num_ant_x, dist_x; c0 = 1)
    positions = manifold.scaled_antenna_positions
    @test get_num_ants(manifold) == num_ant_x^2
    @test get_steer_vec(manifold, [0, 0, 1]) ≈ ones(num_ant_x^2)
    @test get_steer_vec(manifold, [1, 0, 0]) ≈ sign.(positions[1,:]).*im
    @test get_steer_vec(manifold, [0, 1, 0]) ≈ sign.(positions[2,:]).*im
    @test all(abs.(rem.(positions[1,:] .- π/2, π, RoundNearest)) .≈ 0)
    @test all(abs.(rem.(positions[2,:] .- π/2, π, RoundNearest)) .≈ 0)

    num_ant_y = 3
    dist_y = 1
    manifold = IdealManifold(1, num_ant_x, num_ant_y, dist_x, dist_y; c0 = 1)
    positions = manifold.scaled_antenna_positions
    @test get_num_ants(manifold) == num_ant_x * num_ant_y
    @test get_steer_vec(manifold, [0, 0, 1]) ≈ ones(num_ant_x*num_ant_y)
    @test get_steer_vec(manifold, [1, 0, 0]) ≈ sign.(positions[1,:]).*im
    @test get_steer_vec(manifold, [0, 1, 0]) ≈ ones(num_ant_x*num_ant_y)
    @test all(rem.(positions[1,:] .- π/2, π, RoundNearest) .≈ 0)
    @test all(rem.(positions[2,:], 2π, RoundNearest) .≈ 0)
end
