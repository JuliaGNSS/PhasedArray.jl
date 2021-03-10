@testset "Normalize LUT" begin
    lut = SVector(complex.(ones(4,4), ones(4,4)) * 3, complex.(ones(4,4), ones(4,4)) * 2)
    lut_array = [lut[i][j,k] for i = 1:length(lut), j = 1:size(lut[1],1), k = 1:size(lut[1],2)]
    normalized_lut = @inferred PhasedArray.norm_manifold(lut_array)
    lut_test = SVector(complex.(ones(4,4), ones(4,4)) * 3 / sqrt(26 / 2), complex.(ones(4,4), ones(4,4)) * 2 / sqrt(26 / 2))
    lut_array_test = [lut_test[i][j,k] for i = 1:length(lut_test), j = 1:size(lut_test[1],1), k = 1:size(lut_test[1],2)]
    @test normalized_lut ≈ lut_array_test
end

@testset "LUT expansion" begin
    lut = SVector(randn(ComplexF64, 10, 10), randn(ComplexF64, 10, 10), randn(ComplexF64, 10, 10), randn(ComplexF64, 10, 10))
    lut_array = [lut[i][j,k] for i = 1:length(lut), j = 1:size(lut[1],1), k = 1:size(lut[1],2)]
    lut_expanded = @inferred PhasedArray.expand(lut_array, 2)
    num_ants = length(lut)
    @test size(lut_expanded) == (num_ants, 14, 14)
    for i = 1:num_ants
        @test lut_expanded[i,:,1] == lut_expanded[i,:,11]
        @test lut_expanded[i,:,2] == lut_expanded[i,:,12]
        @test lut_expanded[i,:,3] == lut_expanded[i,:,13]
        @test lut_expanded[i,:,4] == lut_expanded[i,:,14]
        @test lut_expanded[i,1,3:12] == circshift(lut_expanded[i,5,3:12], 5)
        @test lut_expanded[i,2,3:12] == circshift(lut_expanded[i,4,3:12], 5)
        @test lut_expanded[i,14,3:12] == circshift(lut_expanded[i,10,3:12], 5)
        @test lut_expanded[i,13,3:12] == circshift(lut_expanded[i,11,3:12], 5)
    end

    @test PhasedArray.calc_expansion_length(Constant) == 0
    @test PhasedArray.calc_expansion_length(Constant) == 0
    @test PhasedArray.calc_expansion_length(Quadratic) == 18
end

@testset "LUT manifold" begin
    λ = 0.1904
    ant_positions = SVector(λ / 4 * SVector(1, 1, 0), λ / 4 * SVector(-1, 1, 0), λ / 4 * SVector(1, -1, 0), λ / 4 * SVector(-1, -1 , 0))
    num_ants = length(ant_positions)
    # Simulated measured pattern
    antenna_gain = 2
    sph2cart = CartesianFromSpherical()
    lut = map(ant_position -> [antenna_gain * cis(2 * π / λ * sph2cart(Spherical(1.0, θ * π / 180, π / 2 - ϕ * π / 180))' * ant_position) for ϕ=0:90, θ=0:359], ant_positions)
    manifold = @inferred RealManifold(lut, max_elevation = π / 2)
    @test typeof(manifold) <: AbstractManifold{4}
    @test norm(@inferred(get_steer_vec(manifold, SVector(1, 1, 1)))) ≈ sqrt(num_ants)
    @test @inferred(get_steer_vec(manifold, SVector(0,0,1))) ≈ [1, 1, 1, 1]
    @test @inferred(get_steer_vec(manifold, SVector(0,1,0))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0))) ≈ [1im, -1im, 1im, -1im]

    @test @inferred(get_steer_vec(manifold, Spherical(1,0,0))) ≈ [1im, -1im, 1im, -1im]
    @test @inferred(get_steer_vec(manifold, Spherical(1,0,0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]

    manifold = @inferred RealManifold(lut[1], lut[2], lut[3], lut[4], max_elevation = π / 2)
    @test typeof(manifold) <: AbstractManifold{4}
    @test norm(@inferred(get_steer_vec(manifold, SVector(1, 1, 1)))) ≈ sqrt(num_ants)
    @test @inferred(get_steer_vec(manifold, SVector(0,0,1))) ≈ [1, 1, 1, 1]
    @test @inferred(get_steer_vec(manifold, SVector(0,1,0))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0))) ≈ [1im, -1im, 1im, -1im]

    @test @inferred(get_steer_vec(manifold, Spherical(1,0,0))) ≈ [1im, -1im, 1im, -1im]
    @test @inferred(get_steer_vec(manifold, Spherical(1,0,0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0), RotXYZ(0.0,0.0,π/2))) ≈ [1im, 1im, -1im, -1im]
end

@testset "LUT manifold interpolated" begin
    λ = 0.1904
    ant_positions = SVector(λ / 4 * SVector(1, 1, 0), λ / 4 * SVector(-1, 1, 0), λ / 4 * SVector(1, -1, 0), λ / 4 * SVector(-1, -1 , 0))
    num_ants = length(ant_positions)
    # Simulated measured pattern
    antenna_gain = 2;
    sph2cart = CartesianFromSpherical()
    lut = map(ant_position -> [antenna_gain * cis(2 * π / λ * sph2cart(Spherical(1.0, θ * π / 180, π / 2 - ϕ * π / 180))' * ant_position) for ϕ=0:90, θ=0:359], ant_positions)
    manifold = @inferred RealManifold(lut, Linear, max_elevation = π / 2)
    @test typeof(manifold) <: AbstractManifold{4}

    @test isapprox(norm(@inferred(get_steer_vec(manifold, SVector(1, 1, 1)))), sqrt(num_ants); rtol = 1e-4)
    @test @inferred(get_steer_vec(manifold, SVector(0,0,1))) ≈ [1, 1, 1, 1]
    @test @inferred(get_steer_vec(manifold, SVector(0,1,0))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(get_steer_vec(manifold, SVector(1,0,0))) ≈ [1im, -1im, 1im, -1im]
end

@testset "Example LUT manifold" begin
    num_ants = length(EXAMPLE_LUT)
    manifold = RealManifold(EXAMPLE_LUT)
    @test typeof(manifold) <: AbstractManifold{2}
    @test norm(@inferred(get_steer_vec(manifold, SVector(0,0,1)))) < sqrt(num_ants)
end