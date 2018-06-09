module PhasedArray

  using Interpolations, CoordinateTransformations
  import Base.transpose

  export manifold, pattern_plotting_data, plot_pattern, pattern_3D_plotting_data, plot_pattern_3D, plot_manifold_3D

  include("manifold.jl")
  include("plots.jl")
end
