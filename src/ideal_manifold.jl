struct IdealManifold{N} <: AbstractManifold{N}
    scaled_antenna_positions::SVector{N,SVector{3,Float64}}
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
function IdealManifold(f_0, antenna_positions::SVector{N,SVector{3,T}}; c0 = 299_792_458) where {N, T <: Real}
    λ = c0 / f_0
    IdealManifold(map(antenna_position -> 2π / λ * antenna_position, antenna_positions))
end

function IdealManifold(f_0, antenna_positions::SVector{3,T}...; c0 = 299_792_458) where T <: Real
    IdealManifold(f_0, SVector(antenna_positions), c0 = c0)
end

function get_steer_vec(manifold::IdealManifold, doa)
    map(scaled_antenna_position -> cis(scaled_antenna_position' * doa), manifold.scaled_antenna_positions)
end

function IdealManifold()
    IdealManifold(0.0, SVector{1, SVector{3, Float64}}(SVector(0.0,0.0,0.0)))
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
