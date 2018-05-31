
"""
  manifold(ant_pos::Array{Float64, 2}, f_0)

Calculates the manifold based on antenna position and signal frequency.
Returns a function to get the steering vector based on DOA.

# Examples
```julia-repl
julia> steer_vec = manifold(0.1904 / 4 * [1 -1 1 -1; 1 1 -1 -1; 0 0 0 0], 1575420e3)
julia> cur_steer_vec = steer_vec(SVector(0,0,1))
```
"""
function manifold(ant_pos::Array{Float64, 2}, f_0, ampl = doa -> 1, c₀ = 299_792_458)
  λ = c₀ / f_0
  doa -> cis.(2 * π / λ * doa' * ant_pos).' .* ampl(doa);
end

"""
  manifold(lut::Array{Complex{Float64}, 3})
  manifold(lut::Array{Complex{Float64}, 3}, interpolation)

Calculates the manifold based on a Look up Table (LuT).
Returns a function to get the steering vector based on DOA.
Interpolation is done by B-splines.
The parameter `interpolation` can be set to Constant(), Linear(), Quadratic(Reflect()), Cubic(Reflect())
For more information see https://github.com/JuliaMath/Interpolations.jl

# Examples
```julia-repl
julia> steer_vec = manifold(lut)
julia> cur_steer_vec = steer_vec(SVector(0,0,1))
```
"""

function manifold(lut::Array{Complex{Float64}, 3}, interpolation = Constant(), max_el = π)
  test_lut_correctness(lut, max_el)
  num_ants = size(lut, 1)
  num_ϕs = size(lut, 2)
  num_θs = size(lut, 3)
  res_ϕ = max_el / (num_ϕs - 1)
  res_θ = 2 * π / num_θs
  lut = norm_manifold(lut)
  lut_expanded = Array{Complex{Float64}, 3}(num_ants, num_ϕs, num_θs + 2)
  lut_expanded[:,:,2:num_θs + 1] = lut
  lut_expanded[:,:,1] = lut[:,:,num_θs]
  lut_expanded[:,:,num_θs + 2] = lut[:,:,1]
  itp = interpolate(lut_expanded, (NoInterp(), BSpline(interpolation), BSpline(interpolation)), OnGrid())
  doa -> _get_steer_vec(doa, itp, res_ϕ, res_θ, num_θs, num_ants, max_el)
end

function _get_steer_vec(doa::Spherical, itp_lut, res_ϕ, res_θ, num_θs, num_ants, max_el)
  ϕ = π / 2 - doa.ϕ # convert to mathematic
  θ = doa.θ + π * ((ϕ < 0) || (ϕ > π))
  ϕ = ϕ - (2 * ϕ) * (ϕ < 0) - (2 * (ϕ - π)) * (ϕ > π && max_el == π) - (ϕ - max_el) * (ϕ > max_el) # 0 <= ϕ <= max_el
  idx_ϕ = ϕ / res_ϕ + 1
  idx_θ = mod(θ / res_θ, num_θs) + 2
  [itp_lut[ant,idx_ϕ,idx_θ]::Complex{Float64} for ant = 1:num_ants]
end

function _get_steer_vec(doa, itp_lut, res_ϕ, res_θ, num_θs, num_ants, max_el)
  _get_steer_vec(Spherical(doa), itp_lut, res_ϕ, res_θ, num_θs, num_ants, max_el)
end

"""
  norm_manifold(lut::Array{Complex{Float64}, 3})

Normalizes the manifold such that the maximal norm is 2.
"""

function norm_manifold(lut::Array{Complex{Float64}, 3})
  num_ants = size(lut, 1)
  steer_vecs_powers = sum(abs2.(lut), 1);
  max_gain = maximum(vec(steer_vecs_powers)) / num_ants;
  lut ./ sqrt(max_gain);
end

function test_lut_correctness(lut::Array{Complex{Float64}, 3}, max_el)
  num_els = size(lut, 2)
  near_horizont_index = num_els * (max_el < π / 2) + floor(Int, π / 2 / max_el * num_els) * (max_el >= π / 2)
  zenith_angle_std = mean(std(angle.(lut[:,1,:]), 1))
  horizont_angle_std = mean(std(angle.(lut[:,near_horizont_index,:]), 1))
  (zenith_angle_std > horizont_angle_std) && error("First row of manifold LUT should be zenith")
  error_az = norm(lut[1,:,1] - lut[1,:,end]) / num_els
  (error_az > 2) && error("Azimuth of manifold LUT should vary over columns. First and last column should be similar.")
end

Base.transpose(sph::Spherical) = CartesianFromSpherical()(sph).'
