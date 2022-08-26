struct RealManifold{N, T <: AbstractInterpolation} <: AbstractManifold{N}
    num_azimuth_angles::Int
    num_elevation_angles::Int
    azimuth_step::Float64
    elevation_step::Float64
    max_elevation::Float64
    expansion_length::Int
    lut::T
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
Possible values are: `Constant`, `Linear` or `Quadratic`. For more information
see regarding the interpolation, see:
https://github.com/JuliaMath/Interpolations.jl.
"""
function RealManifold(
    lut::SVector{N, <: AbstractMatrix{Complex{T}}},
    interpolation::Type{Q} = Constant;
    max_elevation = 1π,
    normalize = true
) where {N, Q <: Interpolations.Degree, T <: Real}
    lut_array = [lut[i][j,k] for i = 1:length(lut), j = 1:size(lut[1],1), k = 1:size(lut[1],2)]
    RealManifold(
        lut_array,
        NumAnts(N),
        Q,
        max_elevation = max_elevation,
        normalize = normalize
    )
end

function RealManifold(
    lut::AbstractArray{Complex{T}, 3},
    num_ants::NumAnts{N},
    interpolation::Type{Q} = Constant;
    max_elevation = 1π,
    normalize = true
) where {N, Q <: Interpolations.Degree, T <: Real}
    test_lut_correctness(lut, max_elevation)
    num_elevation_angles = size(lut, 2)
    num_azimuth_angles = size(lut, 3)
    elevation_step = max_elevation / (num_elevation_angles - 1)
    azimuth_step = 2 * π / num_azimuth_angles
    normalized_lut = normalize ? norm_manifold(lut) : lut
    expansion_length = calc_expansion_length(interpolation)
    expanded_lut = expand(normalized_lut, expansion_length)
    interpolated_lut = interpolate(expanded_lut, (NoInterp(), BSpline(get_elevation_interpolation(Q)), BSpline(get_azimuth_interpolation(Q))))
    RealManifold{N, typeof(interpolated_lut)}(num_azimuth_angles, num_elevation_angles, azimuth_step, elevation_step, max_elevation, expansion_length, interpolated_lut)
end

get_elevation_interpolation(::Type{Q}) where Q <: Interpolations.DegreeBC{0} = Q()
get_elevation_interpolation(::Type{Q}) where Q <: Interpolations.DegreeBC{1} = Q()
# In Elevation direction the boundary is reflective with a 180° circleshift, see `expand`
get_elevation_interpolation(::Type{Q}) where Q <: Interpolations.DegreeBC = Q(Reflect(OnCell()))
get_azimuth_interpolation(::Type{Q}) where Q <: Interpolations.DegreeBC{0} = Q()
get_azimuth_interpolation(::Type{Q}) where Q <: Interpolations.DegreeBC{1} = Q()
# In Azimuth direction the boundary is periodic so one could use Quadratic(Periodic(OnCell())), but 
# this is twice as slow as Quadratic(Reflect(OnCell())), that's why we copy left to the right and right to left
get_azimuth_interpolation(::Type{Q}) where Q <: Interpolations.DegreeBC = Q(Reflect(OnCell()))

function RealManifold(lut::AbstractMatrix{Complex{T}}...; interpolation::Type{Q} = Constant, max_elevation = 1π, normalize = true) where {Q <: Interpolations.Degree, T <: Real}
    RealManifold(SVector(lut), Q, max_elevation = max_elevation, normalize = normalize)
end

function get_steer_vec(manifold::RealManifold{N}, doa::Spherical) where N
    elevation = π / 2 - doa.ϕ # convert to mathematic
    # Rotate azimuth by 180° if elevation is less than zero or greater than π
    azimuth = doa.θ + π * ((elevation < 0) || (elevation > π))
    # Convey the elevation in between 0 <= ϕ <= π
    elevation -= (elevation < 0) * (2 * elevation) +
        (elevation > π && manifold.max_elevation ≈ π) * (2 * (elevation - π))
    elevation > manifold.max_elevation && error("Elevation is above max elevation. Try to increase max elevation.")
    elevation_index = elevation / manifold.elevation_step + 1 + manifold.expansion_length
    azimuth_index = mod(azimuth / manifold.azimuth_step, manifold.num_azimuth_angles) + 1 + manifold.expansion_length
    get_steer_vec(manifold, azimuth_index, elevation_index)
end

function get_steer_vec(manifold::RealManifold{N}, azimuth_index::Real, elevation_index::Real) where N
    steer_vec = MVector{N, ComplexF64}(undef)
    @inbounds for i = 1:length(steer_vec)
        steer_vec[i] = manifold.lut(i, elevation_index, azimuth_index)
    end
    SVector(steer_vec)
end

# See this issue: https://github.com/JuliaGeometry/CoordinateTransformations.jl/issues/25
# if you'd like to know the direction of θ and ϕ
function get_steer_vec(manifold::RealManifold, doa)
    doa_sph = cart2sph(doa) # θ can be between -π and π
    elevation = π / 2 - doa_sph.ϕ # convert to mathematic
    elevation > manifold.max_elevation && error("Elevation is above max elevation. Try to increase max elevation.")
    azimuth = doa_sph.θ + 2π * (doa_sph.θ < 0)
    elevation_index = elevation / manifold.elevation_step + 1 + manifold.expansion_length
    azimuth_index = azimuth / manifold.azimuth_step + 1 + manifold.expansion_length
    get_steer_vec(manifold, azimuth_index, elevation_index)
end

function calc_expansion_length(::Type{<:Interpolations.Degree})
    1
end

function calc_expansion_length(::Type{Quadratic})
    18
end

# Cubic has greater errors than Quadratic when compared with BSpline(Quadratic(Periodic(OnCell())))? Force error
function calc_expansion_length(::Type{Cubic})
    error("Could not find sufficient expansion length. Maybe implement boundary to be periodic?") # 1
end

function expand(lut::AbstractArray{Complex{T}, 3}, num_expand::Int) where {T <: Real}
    num_expand >= 0 || error("Expand number must be greater than zero")
    num_ants = size(lut, 1)
    num_ϕs = size(lut, 2)
    num_θs = size(lut, 3)
    lut_expanded = Array{Complex{T}, 3}(undef, num_ants, num_ϕs + 2 * num_expand, num_θs + 2 * num_expand)
    lut_expanded[:,num_expand + 1:num_ϕs + num_expand,num_expand + 1:num_θs + num_expand] .= lut
    # In Elevation direction the boundary is reflective with a 180° circleshift
    for ant = 1:num_ants
        for i = 1:num_expand
            lut_expanded[ant,num_expand - i + 1,num_expand + 1:num_θs + num_expand] .= circshift(lut_expanded[ant,num_expand + i + 1,num_expand + 1:num_θs + num_expand], floor(Int, num_θs / 2))
            lut_expanded[ant,num_ϕs + num_expand + i,num_expand + 1:num_θs + num_expand] .= circshift(lut_expanded[ant,num_ϕs + num_expand - i,num_expand + 1:num_θs + num_expand], floor(Int, num_θs / 2))
        end
        # In Azimuth direction the boundary is periodic so one could use Quadratic(Periodic(OnCell())), but 
        # this is twice as slow as Quadratic(Reflect(OnCell())), that's why we copy left to the right and right to left
        for i = 1:num_expand
            lut_expanded[ant,:,num_θs + num_expand + i] .= lut_expanded[ant,:,num_expand + i]
            lut_expanded[ant,:,num_expand - i + 1] .= lut_expanded[ant,:,num_θs + num_expand - i + 1]
        end
    end
    lut_expanded
end

"""
  norm_manifold(lut::Array{Complex{Float64}, 3})

Normalizes the manifold such that the maximal norm is `sqrt(num_ants)`.
"""
function norm_manifold(lut::AbstractArray{Complex{T}, 3}) where {T <: Real}
    max_norm = mapreduce(norm, max, Slices(lut, True(), False(), False()))
    lut ./ (max_norm / sqrt(size(lut, 1)))
end

function test_lut_correctness(lut::AbstractArray{Complex{T}, 3}, max_el) where T <: Real
    num_els = size(lut, 2)
    num_ants = size(lut, 1)
    near_horizont_index = min(floor(Int, π / 2 / max_el * num_els), num_els)
    lut_reference_one = lut[2:end, :, :] ./ lut[1:1, :, :]
    zenith_angle_var = mean(var(angle.(@view(lut_reference_one[:, 1, :])), dims = 2))
    horizont_angle_var = mean(var(angle.(@view(lut_reference_one[:, near_horizont_index, :])), dims = 2))
    (zenith_angle_var > horizont_angle_var) && error("First row of manifold LUT should be zenith")
    error_az = norm(lut[:,:,1] - lut[:,:,end]) / num_els / num_ants
    (error_az > 2) && error("Azimuth of manifold LUT should vary over columns. First and last column should be similar.")
end
