@testset "DOA MUSIC" begin
    ant_pos = @SMatrix([-1 1 -1 1; -1 -1 1 1; 0 0 0 0])
    manifold = IdealManifold(1, 1 / 4 * ant_pos; c0 = 1)

    @test doa_music(manifold, SVector( 1,  1,  1,  1)) == Spherical(1.0, 0, π/2)
    @test doa_music(manifold, SVector( 1,  1, -1, -1)) == Spherical(1.0, π/2, 0)
    @test doa_music(manifold, SVector(-1,  1, -1,  1)) == Spherical(1.0, 0, 0)

    doa = Spherical(1, 137π/180, 40π/180)
    @test doa_music(manifold, get_steer_vec(manifold, doa)) == doa
end
