module PhasedArray

    using Interpolations, CoordinateTransformations, LinearAlgebra, StaticArrays, Statistics, PGFPlotsX
    import Base.transpose, Base.filter

    export
        manifold,
        Pattern,
        Pattern3D,
        plot,
        filter,
        calc_prewhitening_filter,
        calc_amplitude_filter,
        calc_eigen_beamformer,
        calc_variance_covariance

    const cart2sph = SphericalFromCartesian()
    const sph2cart = CartesianFromSpherical()

    include("manifold.jl")
    include("pattern.jl")
    include("pgfplots.jl")
    include("filter.jl")
end
