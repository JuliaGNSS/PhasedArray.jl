@testset "SAMV" begin
    f_0 = 1575420e3
    c₀ = 299_792_458
    λ = c₀ / f_0
    ant_pos = SVector(λ / 4 * SVector(1, 1, 0), λ / 4 * SVector(-1, 1, 0), λ / 4 * SVector(1, -1, 0), λ / 4 * SVector(-1, -1 , 0))
    num_ants = length(ant_pos)
    manifold = IdealManifold(f_0, ant_pos)

    signal1_doa_sph = Spherical(1.0, 0.1π, 0.8π / 2)
    signal2_doa_sph = Spherical(1.0, 1.2π, 0.5π / 2)

    signal1_doa = CartesianFromSpherical()(signal1_doa_sph)
    signal2_doa = CartesianFromSpherical()(signal2_doa_sph)

    num_samples = 1000
    signal1_freq = 1.2e6
    signal2_freq = -0.5e6
    sampling_freq = 5e6
    signal = transpose(get_steer_vec(manifold, signal1_doa_sph)) .* cis.(2π * (0:num_samples - 1) * signal1_freq / sampling_freq .+ rand() * 2π) * 15.0 +
        transpose(get_steer_vec(manifold, signal2_doa_sph)) .* cis.(2π * (0:num_samples - 1) * signal2_freq / sampling_freq .+ rand() * 2π) * 20.0 +
        randn(ComplexF64, num_samples, num_ants)

    Rxx = Hermitian(transpose(signal) * conj(signal) / num_samples)

    doas, signal_powers = PhasedArray.samv2(manifold, Rxx; num_doas = 1000,)

    found_doas = doas[signal_powers .> 13^2]
    found_signal_powers = signal_powers[signal_powers .> 13^2]
    sort_perm = sortperm(found_signal_powers)

    @test acosd(signal2_doa' * found_doas[sort_perm][2]) < 1 # Smaller than 1 degree
    @test acosd(signal1_doa' * found_doas[sort_perm][1]) < 1.5 # Smaller than 1.5 degree
    
#=
    pattern = ScatterPattern(doas, sqrt.(abs.(signal_powers)))
    plot(pattern)
    scatter!([signal1_doa.θ], [(π / 2 - signal1_doa.ϕ) * 180 / π], marker = :cross, markerstrokewidth = 4)
    scatter!([signal2_doa.θ], [(π / 2 - signal2_doa.ϕ) * 180 / π], marker = :cross, markerstrokewidth = 4)
=#
end