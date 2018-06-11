module PhasedArray

  using Interpolations, CoordinateTransformations, PyCall, PyPlot
  import Base.transpose

  export manifold, pattern_plotting_data, plot_pattern, pattern_3D_plotting_data, plot_pattern_3D, plot_manifold_3D, draw_polar_axes, init_animate_pattern_data, animate_pattern

  include("manifold.jl")
  include("plots.jl")
  include("animate_pattern.jl")
end
