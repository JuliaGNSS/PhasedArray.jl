struct RealManifold{N, T <: AbstractInterpolation} <: AbstractManifold{N}
    num_azimuth_angles::Int
    num_elevation_angles::Int
    azimuth_step::Float64
    elevation_step::Float64
    max_elevation::Float64
    expansion_length::Int
    lut::SVector{N, T}
end

"""
LUT must be a vector with the length of the number of antenna elements of type
SVector of Matrices.
Each matrix must have the size of #elevation angles x #azimuth angles. The
azimuth angle starts from the x-Axis and increases towards the y-Axis. The
zenith must be at the top of the matrix and the horizont / bottom must be at
the bottom of the matrix. The maximal elevation is given by `max_elevation`.
A value of `π/2` means horizont and `π` until bottom of antenna.
The type of interpolation can be set by the optional parameter `interpolation`.
Possible values are: `Constant()`, `Linear()`, `Quadratic(Reflect(OnCell()))`
or `Cubic(Reflect(OnCell()))`. For more information see regarding the
interpolation, see: https://github.com/JuliaMath/Interpolations.jl.
"""
function RealManifold(lut::SVector{N, Array{Complex{Float64}, 2}}, interpolation = Constant(), max_elevation = 1π) where N
    test_lut_correctness(lut, max_elevation)
    num_elevation_angles = size(lut[1], 1)
    num_azimuth_angles = size(lut[1], 2)
    elevation_step = max_elevation / (num_elevation_angles - 1)
    azimuth_step = 2 * π / num_azimuth_angles
    normalized_lut = norm_manifold(lut)
    expansion_length = calc_expansion_length(interpolation)
    expanded_lut = expand(normalized_lut, expansion_length)
    interpolated_lut = map(X -> interpolate(X, (BSpline(interpolation), BSpline(interpolation))), expanded_lut)
    RealManifold(num_azimuth_angles, num_elevation_angles, azimuth_step, elevation_step, max_elevation, expansion_length, interpolated_lut)
end

function get_steer_vec(manifold::RealManifold, doa::Spherical)
    elevation = π / 2 - doa.ϕ # convert to mathematic
    # Rotate azimuth by 180° if elevation is less than zero or greater than π
    azimuth = doa.θ + π * ((elevation < 0) || (elevation > π))
    # Convey the elevation in between 0 <= ϕ <= max_elevation
    elevation -= (elevation < 0) * (2 * elevation)
        + (elevation > π && manifold.max_elevation ≈ π) * (2 * (elevation - π))
        + (elevation > manifold.max_elevation) * (elevation - manifold.max_elevation)
    elevation_index = elevation / manifold.elevation_step + 1 + manifold.expansion_length
    azimuth_index = mod(azimuth / manifold.azimuth_step, manifold.num_azimuth_angles) + 1 + manifold.expansion_length
    map(X -> X(elevation_index,azimuth_index), manifold.lut)
end

function get_steer_vec(manifold::RealManifold, doa)
    get_steer_vec(manifold, cart2sph(doa))
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

function expand(lut::SVector{N, Array{Complex{Float64}, 2}}, num_expand::Int) where N
    num_expand >= 0 || error("Expand number must be greater than zero")
    num_ϕs = size(lut[1], 1)
    num_θs = size(lut[1], 2)
    lut_expanded = map(X -> Matrix{eltype(eltype(lut))}(undef, num_ϕs + 2 * num_expand, num_θs + 2 * num_expand), lut)
    foreach((X, Y) -> X[num_expand + 1:num_ϕs + num_expand,num_expand + 1:num_θs + num_expand] .= Y, lut_expanded, lut)
    for ant = 1:N
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

"""
  norm_manifold(lut::Array{Complex{Float64}, 3})

Normalizes the manifold such that the maximal norm is `sqrt(num_ants)`.
"""
function norm_manifold(lut::SVector{N, Array{Complex{Float64}, 2}}) where N
    max_gain = 0.0
    length_lut = length(lut[1])
    for j = 1:length_lut
        current_gain = 0.0
        for i = 1:N
            current_gain += abs2(lut[i][j])
        end
        max_gain = max(max_gain, current_gain)
    end
    map(X -> X ./ sqrt(max_gain / N), lut)
end

function test_lut_correctness(lut::SVector{N, Array{Complex{Float64}, 2}}, max_el) where N
    num_els = size(lut[1], 1)
    near_horizont_index = num_els * (max_el < π / 2) + floor(Int, π / 2 / max_el * num_els) * (max_el >= π / 2)
    lut_reference_one = map(X -> X ./ lut[1], lut[2:end])
    zenith_angle_var = mapreduce(X -> var(angle.(X[1, :])), +, lut_reference_one)
    horizont_angle_var = mapreduce(X -> var(angle.(X[near_horizont_index, :])), +, lut_reference_one)
    (zenith_angle_var > horizont_angle_var) && error("First row of manifold LUT should be zenith")
    error_az = norm(lut[1][:,1] - lut[1][:,end]) / num_els
    (error_az > 2) && error("Azimuth of manifold LUT should vary over columns. First and last column should be similar.")
end
