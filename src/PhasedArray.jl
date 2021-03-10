module PhasedArray

    using
        Interpolations,
        CoordinateTransformations,
        LinearAlgebra,
        StaticArrays,
        Statistics,
        JuliennedArrays

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
        calc_variance_covariance

    const cart2sph = SphericalFromCartesian()
    const sph2cart = CartesianFromSpherical()

    struct NumAnts{N} end
    NumAnts(N) = NumAnts{N}()

    include("manifold.jl")
    include("ideal_manifold.jl")
    include("real_manifold.jl")
    include("pattern.jl")
    include("filter.jl")
end
