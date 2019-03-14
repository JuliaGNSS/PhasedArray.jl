module PhasedArray

    using Interpolations, CoordinateTransformations, Statistics, LinearAlgebra, Printf, PyCall, PyPlot
    import Base.transpose, Base.filter

    export
        manifold,
        pattern_plotting_data,
        plot_pattern,
        pattern_3D_plotting_data,
        plot_pattern_3D,
        plot_manifold_3D,
        draw_polar_axes,
        init_animate_pattern_data,
        animate_pattern,
        filter,
        calc_prewhitening_filter,
        calc_amplitude_filter,
        calc_eigen_beamformer,
        calc_variance_covariance

    include("manifold.jl")
    include("plots.jl")
    #include("animate_pattern.jl")
    include("filter.jl")
end
