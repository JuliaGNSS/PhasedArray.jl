module PhasedArray

    using Interpolations, CoordinateTransformations, LinearAlgebra, StaticArrays, Statistics, PGFPlotsX
    import Base.transpose, Base.filter

    export
        IdealManifold,
        RealManifold,
        get_steer_vec,
        Pattern,
        Pattern3D,
        filter,
        calc_prewhitening_filter,
        calc_amplitude_filter,
        calc_eigen_beamformer,
        calc_variance_covariance

    const cart2sph = SphericalFromCartesian()
    const sph2cart = CartesianFromSpherical()

    include("manifold.jl")
    include("ideal_manifold.jl")
    include("real_manifold.jl")
    include("pattern.jl")
    include("pgfplots.jl")
    include("filter.jl")
end
