@testset "Pattern plotting data" begin
    f_0 = 1575420e3
    c₀ = 299_792_458
    λ = c₀ / f_0
    ant_pos = λ / 4 * [1 -1 1 -1; 1 1 -1 -1; 0 0 0 0];
    get_steer_vec = PhasedArray.manifold(ant_pos, f_0)

    azs, els, values = @inferred pattern_plotting_data(get_steer_vec, norm, 360, 91)
    @test length(azs) == 360
    @test length(els) == 91
    @test size(values) == (91, 360)
end

@testset "3D Pattern plotting data" begin
    f_0 = 1575420e3
    c₀ = 299_792_458
    λ = c₀ / f_0
    ant_pos = λ / 4 * [1 -1 1 -1; 1 1 -1 -1; 0 0 0 0];
    get_steer_vec = PhasedArray.manifold(ant_pos, f_0)

    X, Y, Z, gains = @inferred pattern_3D_plotting_data(get_steer_vec, norm, 360, 181, 180)
    @test size(X) == (181, 360)
    @test size(Y) == (181, 360)
    @test size(Z) == (181, 360)
    @test size(gains) == (181, 360)
end