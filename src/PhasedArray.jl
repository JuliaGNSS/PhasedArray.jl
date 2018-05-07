module PhasedArray

  using Interpolations, CoordinateTransformations, PyPlot
  import Base.transpose

  export manifold, plot_pattern, plot_pattern_3D, plot_manifold_3D

  include("manifold.jl")
  include("plots.jl")
end
