
"""
  manifold(ant_pos::Array{Float64, 2}, f_0)

Calculates the manifold based on antenna position and signal frequency.
Returns a function to get the steering vector based on DOA.

# Examples
```julia-repl
julia> steer_vec = manifold(0.1904 / 4 * SVector(SVector(1, 1, 0), SVector(-1, 1, 0), SVector(1, -1, 0), SVector(-1, -1, 0)), 1575420e3)
julia> cur_steer_vec = steer_vec(SVector(0,0,1))
```
"""
function manifold(antenna_positions::SVector{Q,SVector{3,T}}, f_0, c₀ = 299_792_458) where Q where T <: Real
    λ = c₀ / f_0
    scaled_antenna_positions = map(antenna_position -> 2π / λ * antenna_position, antenna_positions)
    doa -> _get_steer_vec(doa, scaled_antenna_positions)
end

function _get_steer_vec(doa, scaled_antenna_positions::SVector{Q,SVector{3,Float64}}) where Q
    map(scaled_antenna_position -> cis(scaled_antenna_position' * doa), scaled_antenna_positions)
end

function _get_steer_vec(doa::Spherical, scaled_antenna_positions::SVector{Q,SVector{3,Float64}}) where Q
    _get_steer_vec(sph2cart(doa), scaled_antenna_positions)
end

"""
  manifold(lut::SVector{Q, Array{Complex{Float64}, 2}}) where Q
  manifold(lut::SVector{Q, Array{Complex{Float64}, 2}}, interpolation) where Q

Calculates the manifold based on a Look up Table (LuT).
Returns a function to get the steering vector based on DOA.
Interpolation is done by B-splines.
The parameter `interpolation` can be set to Constant(), Linear(), Quadratic(Reflect(OnCell())), Cubic(Reflect(OnCell()))
For more information see https://github.com/JuliaMath/Interpolations.jl

# Examples
```julia-repl
julia> steer_vec = manifold(lut)
julia> cur_steer_vec = steer_vec(SVector(0,0,1))
```
"""
function manifold(lut::SVector{Q, Array{Complex{Float64}, 2}}, interpolation = Constant(), max_el = π) where Q
    test_lut_correctness(lut, max_el)
    num_ϕs = size(lut[1], 1)
    num_θs = size(lut[1], 2)
    res_ϕ = max_el / (num_ϕs - 1)
    res_θ = 2 * π / num_θs
    lut = norm_manifold(lut)
    expansion_length = calc_expansion_length(interpolation)
    lut_expanded = expand(lut, expansion_length)
    itp = map(X -> interpolate(X, (BSpline(interpolation), BSpline(interpolation))), lut_expanded)
    doa -> _get_steer_vec(doa, itp, res_ϕ, res_θ, num_θs, max_el, expansion_length)
end

function calc_expansion_length(::Constant)
    1
end

function calc_expansion_length(::Linear)
    1
end

function calc_expansion_length(::Quadratic)
    7
end

function expand(lut::SVector{Q, Array{Complex{Float64}, 2}}, num_expand::Int) where Q
    num_expand >= 0 || error("Expand number must be greater than zero")
    num_ϕs = size(lut[1], 1)
    num_θs = size(lut[1], 2)
    lut_expanded = map(X -> Matrix{eltype(eltype(lut))}(undef, num_ϕs + 2 * num_expand, num_θs + 2 * num_expand), lut)
    foreach((X, Y) -> X[num_expand + 1:num_ϕs + num_expand,num_expand + 1:num_θs + num_expand] .= Y, lut_expanded, lut)
    for ant = 1:Q
        for i = 1:num_expand
            lut_expanded[ant][num_expand - i + 1,num_expand + 1:num_θs + num_expand] .= circshift(lut_expanded[ant][num_expand + i + 1,num_expand + 1:num_θs + num_expand], floor(Int, num_θs / 2))
            lut_expanded[ant][num_ϕs + num_expand + i,num_expand + 1:num_θs + num_expand] .= circshift(lut_expanded[ant][num_ϕs + num_expand - i,num_expand + 1:num_θs + num_expand], floor(Int, num_θs / 2))
        end
        for i = 1:num_expand
            lut_expanded[ant][:,num_θs + num_expand + i] .= lut_expanded[ant][:,num_expand + i]
            lut_expanded[ant][:,num_expand - i + 1] .= lut_expanded[ant][:,num_θs + num_expand - i + 1]
        end
    end
    lut_expanded
end

function _get_steer_vec(doa::Spherical, itp_lut, res_ϕ, res_θ, num_θs, max_el, expansion_length)
    ϕ = π / 2 - doa.ϕ # convert to mathematic
    θ = doa.θ + π * ((ϕ < 0) || (ϕ > π))
    ϕ = ϕ - (2 * ϕ) * (ϕ < 0) - (2 * (ϕ - π)) * (ϕ > π && max_el == π) - (ϕ - max_el) * (ϕ > max_el) # 0 <= ϕ <= max_el
    idx_ϕ = ϕ / res_ϕ + 1 + expansion_length
    idx_θ = mod(θ / res_θ, num_θs) + 1 + expansion_length
    map(X -> X(idx_ϕ,idx_θ), itp_lut)
end

function _get_steer_vec(doa, itp_lut, res_ϕ, res_θ, num_θs, max_el, expansion_length)
    _get_steer_vec(cart2sph(doa), itp_lut, res_ϕ, res_θ, num_θs, max_el, expansion_length)
end

"""
  norm_manifold(lut::Array{Complex{Float64}, 3})

Normalizes the manifold such that the maximal norm is `sqrt(num_ants)`.
"""
function norm_manifold(lut::SVector{Q, Array{Complex{Float64}, 2}}) where Q
    max_gain = 0.0
    length_lut = length(lut[1])
    for j = 1:length_lut
        current_gain = 0.0
        for i = 1:Q
            current_gain += abs2(lut[i][j])
        end
        max_gain = max(max_gain, current_gain)
    end
    map(X -> X ./ sqrt(max_gain / Q), lut)
end

function test_lut_correctness(lut::SVector{Q, Array{Complex{Float64}, 2}}, max_el) where Q
    num_els = size(lut[1], 1)
    near_horizont_index = num_els * (max_el < π / 2) + floor(Int, π / 2 / max_el * num_els) * (max_el >= π / 2)
    lut_reference_one = map(X -> X ./ lut[1], lut[2:end])
    zenith_angle_var = mapreduce(X -> var(angle.(X[1, :])), +, lut_reference_one)
    horizont_angle_var = mapreduce(X -> var(angle.(X[near_horizont_index, :])), +, lut_reference_one)
    (zenith_angle_var > horizont_angle_var) && error("First row of manifold LUT should be zenith")
    error_az = norm(lut[1][:,1] - lut[1][:,end]) / num_els
    (error_az > 2) && error("Azimuth of manifold LUT should vary over columns. First and last column should be similar.")
end
