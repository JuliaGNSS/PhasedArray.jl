struct IdealManifold{N, T <: SMatrix{3,N}} <: AbstractManifold{N}
    scaled_antenna_positions::T
end

"""
Calculates the manifold based on antenna position and signal frequency.

# Examples
```julia-repl
julia> manifold = IdealManifold(1575420e3, 0.1904 / 4 * SVector(SVector(1, 1, 0), SVector(-1, 1, 0), SVector(1, -1, 0), SVector(-1, -1, 0)))
julia> get_steer_vec(manifold, SVector(0.0,0.0,1.0))
julia> get_steer_vec(manifold, SVector(0.0,0.0,1.0), RotXYZ(0.0,0.0,0.0))
```
"""
function IdealManifold(f_0::Real, antenna_positions::SMatrix{3,N,T}; c0 = 299_792_458) where {N, T <: Real}
    λ = c0 / f_0
    IdealManifold(2π / λ * antenna_positions)
end

function IdealManifold(f_0::Real, antenna_positions::SVector{N,SVector{3,T}}; c0 = 299_792_458) where {N, T <: Real}
    IdealManifold(f_0, reduce(hcat, antenna_positions), c0 = c0)
end

function IdealManifold(f_0::Real, antenna_positions::SVector{3,T}...; c0 = 299_792_458) where T <: Real
    IdealManifold(f_0, SVector(antenna_positions), c0 = c0)
end

function IdealManifold(f_0::Real, antenna_positions::SVector{3,T}; c0 = 299_792_458) where T <: Real
    IdealManifold(f_0, SMatrix{3,1,T}(antenna_positions), c0 = c0)
end

function IdealManifold(f_0::Real, num_ant_x, num_ant_y, dist_x, dist_y; kwargs...)
    max_x = (num_ant_x - 1) / 2
    max_y = (num_ant_y - 1) / 2
    pos_x = (-max_x:max_x) .* dist_x
    pos_y = (-max_y:max_y) .* dist_y
    ant_pos = SMatrix{3, num_ant_x * num_ant_y}(
        vcat(
            repeat(pos_x, outer = num_ant_y)',
            repeat(pos_y, inner = num_ant_x)',
            zeros(1, num_ant_x * num_ant_y)
        )
    )
    IdealManifold(f_0, ant_pos; kwargs...)
end

function IdealManifold(f_0::Real, num_ant_xy, dist_xy; kwargs...)
    IdealManifold(f_0, num_ant_xy, num_ant_xy, dist_xy, dist_xy; kwargs...)
end

function get_steer_vec(manifold::IdealManifold, doa)
    cis.(transpose(manifold.scaled_antenna_positions) * doa)
end

function IdealManifold()
    IdealManifold(zeros(SMatrix{3, 1, Float64}))
end

function get_steer_vec(manifold::IdealManifold{1}, doa)
    1.0
end

function get_steer_vec(manifold::IdealManifold{1}, doa, attitude)
    1.0
end

function get_steer_vec(manifold::IdealManifold, doa::Spherical)
    get_steer_vec(manifold, sph2cart(doa))
end
