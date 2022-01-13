module PhasedArray

    using
        Interpolations,
        CoordinateTransformations,
        DocStringExtensions,
        LinearAlgebra,
        StaticArrays,
        Statistics,
        JuliennedArrays,
        RecipesBase

    export
        IdealManifold,
        RealManifold,
        AbstractManifold,
        get_steer_vec,
        Pattern,
        Pattern3D,
        calc_prewhitening_filter,
        calc_amplitude_filter,
        calc_eigen_beamformer,
        calc_variance_covariance,
        est_doa,
        est_doa_by_signal_subspace,
        est_doa_by_noise_subspace,
        est_doa_by_music,
        get_num_ants

    const cart2sph = SphericalFromCartesian()
    const sph2cart = CartesianFromSpherical()

    struct NumAnts{N} end
    NumAnts(N) = NumAnts{N}()

    include("manifold.jl")
    include("est_doa.jl")
    include("ideal_manifold.jl")
    include("real_manifold.jl")
    include("pattern.jl")
    include("filter.jl")
end
