struct IdealManifold{Q} <: AbstractManifold
    scaled_antenna_positions::SVector{Q,SVector{3,Float64}}
end

"""
Calculates the manifold based on antenna position and signal frequency.

# Examples
```julia-repl
julia> manifold = IdealManifold(0.1904 / 4 * SVector(SVector(1, 1, 0), SVector(-1, 1, 0), SVector(1, -1, 0), SVector(-1, -1, 0)), 1575420e3)
julia> get_steer_vec(manifold, SVector(0.0,0.0,1.0))
julia> get_steer_vec(manifold, SVector(0.0,0.0,1.0), RotXYZ(0.0,0.0,0.0))
```
"""
function IdealManifold(antenna_positions::SVector{Q,SVector{3,T}}, f_0, c₀ = 299_792_458) where Q where T <: Real
    λ = c₀ / f_0
    IdealManifold(map(antenna_position -> 2π / λ * antenna_position, antenna_positions))
end

function get_steer_vec(manifold::IdealManifold, doa)
    map(scaled_antenna_position -> cis(scaled_antenna_position' * doa), manifold.scaled_antenna_positions)
end

function IdealManifold()
    IdealManifold(SVector{1, SVector{3, Float64}}(SVector(0.0,0.0,0.0)), 0.0)
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
