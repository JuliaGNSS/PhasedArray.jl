@testset "Hill climbing" begin
    values = sin.(0:0.01:2)
    _, idx = @inferred PhasedArray.hill_climbing(values, 10)
    @test idx == CartesianIndex(158)
    values_2D = values * values'
    _, indices = @inferred PhasedArray.hill_climbing(values_2D, 10, 10)
    @test indices == CartesianIndex(158, 158)
end

@testset "Estimate DOA" begin
    ant_pos = @SMatrix([-1 1 -1 1; -1 -1 1 1; 0 0 0 0])
    manifold = IdealManifold(1, 1 / 4 * ant_pos; c0 = 1)

    @test est_doa(manifold, a -> abs(SVector( 1,  1,  1,  1)' * a) / norm(a)) == Spherical(1.0, 0.0, π / 2)
    @test est_doa(manifold, a -> abs(SVector( 1,  1, -1, -1)' * a) / norm(a)) == Spherical(1.0, π / 2, 0.0)
    @test est_doa(manifold, a -> abs(SVector(-1,  1, -1,  1)' * a) / norm(a)) == Spherical(1.0, 0.0, 0.0)

    doa = Spherical(1, 137π / 180, 40π / 180)
    @test est_doa(manifold, a -> abs(get_steer_vec(manifold, doa)' * a) / norm(a)) == doa


    @test PhasedArray.sph2cart(est_doa(manifold, a -> abs(SVector(1, 1, 1, 1)' * a) / norm(a), init_az = 0.1, init_el = π / 2 - 0.1)) ≈ SVector(0, 0, 1)
    @test est_doa(manifold, a -> abs(SVector( 1,  1, -1, -1)' * a) / norm(a), init_az = π / 2 + 0.1, init_el = 0.1) == Spherical(1.0, π / 2, 0.0)
    @test est_doa(manifold, a -> abs(SVector(-1,  1, -1,  1)' * a) / norm(a), init_az = 0.1, init_el = 0.1) == Spherical(1.0, 0.0, 0.0)

    doa = Spherical(1, 137π / 180, 40π / 180)
    @test est_doa(manifold, a -> abs(get_steer_vec(manifold, doa)' * a) / norm(a), init_az = 137π / 180 + 0.1, init_el = 40π / 180 - 0.1) == doa

    @test est_doa_by_signal_subspace(manifold, SVector( 1,  1,  1,  1)) == Spherical(1.0, 0.0, π / 2)
    @test est_doa_by_signal_subspace(manifold, SVector( 1,  1, -1, -1)) == Spherical(1.0, π / 2, 0.0)
    @test est_doa_by_signal_subspace(manifold, SVector(-1,  1, -1,  1)) == Spherical(1.0, 0.0, 0.0)

    doa = Spherical(1, 137π / 180, 40π / 180)
    @test est_doa_by_signal_subspace(manifold, get_steer_vec(manifold, doa)) == doa

    steer_vec = SVector( 1,  1,  1,  1)
    noise_subspace = I - steer_vec * steer_vec' / (steer_vec' * steer_vec)
    @test est_doa_by_noise_subspace(manifold, noise_subspace) == Spherical(1.0, 0.0, π / 2)
    @test est_doa_by_music(manifold, noise_subspace) == Spherical(1.0, 0.0, π / 2)
    steer_vec = SVector( 1,  1, -1, -1)
    noise_subspace = I - steer_vec * steer_vec' / (steer_vec' * steer_vec)
    @test est_doa_by_noise_subspace(manifold, noise_subspace) == Spherical(1.0, π / 2, 0.0)
    @test est_doa_by_music(manifold, noise_subspace) == Spherical(1.0, π / 2, 0.0)

    doa = Spherical(1, 137π / 180, 40π / 180)
    steer_vec = get_steer_vec(manifold, doa)
    noise_subspace = I - steer_vec * steer_vec' / (steer_vec' * steer_vec)
    @test est_doa_by_noise_subspace(manifold, noise_subspace) == doa
    @test est_doa_by_music(manifold, noise_subspace) == doa
end
