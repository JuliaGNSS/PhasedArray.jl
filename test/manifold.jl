@testset "Ideal manifold" begin
    f_0 = 1575420e3
    c₀ = 299_792_458
    λ = c₀ / f_0
    ant_pos = SVector(λ / 4 * SVector(1, 1, 0), λ / 4 * SVector(-1, 1, 0), λ / 4 * SVector(1, -1, 0), λ / 4 * SVector(-1, -1 , 0))
    steer_vec = @inferred manifold(ant_pos, f_0)

    @test @inferred(steer_vec(SVector(0,0,1))) ≈ [1, 1, 1, 1]
    @test @inferred(steer_vec(SVector(0,1,0))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(steer_vec(SVector(1,0,0))) ≈ [1im, -1im, 1im, -1im]
end

@testset "LUT manifold" begin
    λ = 0.1904
    ant_positions = SVector(λ / 4 * SVector(1, 1, 0), λ / 4 * SVector(-1, 1, 0), λ / 4 * SVector(1, -1, 0), λ / 4 * SVector(-1, -1 , 0))
    num_ants = length(ant_positions)
    # Simulated measured pattern
    antenna_gain = 2
    sph2cart = CartesianFromSpherical()
    lut = map(ant_position -> [antenna_gain * cis(2 * π / λ * sph2cart(Spherical(1.0, θ * π / 180, π / 2 - ϕ * π / 180))' * ant_position) for ϕ=0:90, θ=0:359], ant_positions)
    steer_vec = @inferred manifold(lut, Constant(), π / 2)
    @test norm(@inferred(steer_vec(SphericalFromCartesian()(SVector(1, 1, 1))))) ≈ sqrt(num_ants)
    @test @inferred(steer_vec(SphericalFromCartesian()(SVector(0,0,1)))) ≈ [1, 1, 1, 1]
    @test @inferred(steer_vec(SphericalFromCartesian()(SVector(0,1,0)))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(steer_vec(SphericalFromCartesian()(SVector(1,0,0)))) ≈ [1im, -1im, 1im, -1im]
end

@testset "LUT manifold interpolated" begin
    λ = 0.1904
    ant_positions = SVector(λ / 4 * SVector(1, 1, 0), λ / 4 * SVector(-1, 1, 0), λ / 4 * SVector(1, -1, 0), λ / 4 * SVector(-1, -1 , 0))
    num_ants = length(ant_positions)
    # Simulated measured pattern
    antenna_gain = 2;
    sph2cart = CartesianFromSpherical()
    lut = map(ant_position -> [antenna_gain * cis(2 * π / λ * sph2cart(Spherical(1.0, θ * π / 180, π / 2 - ϕ * π / 180))' * ant_position) for ϕ=0:90, θ=0:359], ant_positions)
    steer_vec = @inferred manifold(lut, Linear(), π / 2)

    @test isapprox(norm(@inferred(steer_vec(SphericalFromCartesian()(SVector(1, 1, 1))))), sqrt(num_ants); rtol = 1e-4)
    @test @inferred(steer_vec(SphericalFromCartesian()(SVector(0,0,1)))) ≈ [1, 1, 1, 1]
    @test @inferred(steer_vec(SphericalFromCartesian()(SVector(0,1,0)))) ≈ [1im, 1im, -1im, -1im]
    @test @inferred(steer_vec(SphericalFromCartesian()(SVector(1,0,0)))) ≈ [1im, -1im, 1im, -1im]
end

@testset "Example LUT manifold" begin
    num_ants = length(EXAMPLE_LUT)
    steer_vec = manifold(EXAMPLE_LUT)
    @test norm(@inferred(steer_vec(SphericalFromCartesian()(SVector(0,0,1))))) < sqrt(num_ants)
end

@testset "LUT expansion" begin
    lut = SVector(randn(ComplexF64, 10, 10), randn(ComplexF64, 10, 10), randn(ComplexF64, 10, 10), randn(ComplexF64, 10, 10))
    lut_expanded = PhasedArray.expand(lut, 2)
    num_ants = length(lut)
    for i = 1:num_ants
        @test size(lut_expanded[i]) == (14, 14)
        @test lut_expanded[i][:,1] == lut_expanded[i][:,11]
        @test lut_expanded[i][:,2] == lut_expanded[i][:,12]
        @test lut_expanded[i][:,3] == lut_expanded[i][:,13]
        @test lut_expanded[i][:,4] == lut_expanded[i][:,14]
        @test lut_expanded[i][1,3:12] == circshift(lut_expanded[i][5,3:12], 5)
        @test lut_expanded[i][2,3:12] == circshift(lut_expanded[i][4,3:12], 5)
        @test lut_expanded[i][14,3:12] == circshift(lut_expanded[i][10,3:12], 5)
        @test lut_expanded[i][13,3:12] == circshift(lut_expanded[i][11,3:12], 5)
    end

    @test PhasedArray.calc_expansion_length(Constant()) == 1
    @test PhasedArray.calc_expansion_length(Constant()) == 1
    @test PhasedArray.calc_expansion_length(Quadratic(Reflect(OnCell()))) == 7
end
