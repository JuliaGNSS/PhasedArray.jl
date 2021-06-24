"""
$(SIGNATURES)

Calculate the direction of arrival of the signal `s` to the manifold `m`. The
calculated DOA is returned as `Spherical` object, where the azimuth is counted
counter-clockwise from the x axis and the elevation is counted from the xy-plane
towards the zenith.

# Arguments:
- `m::AbstractManifold{N}`: Manfiold the signal is received with
- `s::AbstractVector`: Signal subspace for which the doa is estimated
"""
function doa_music(m::AbstractManifold{N}, s::AbstractVector; kwargs...) where N
    length(s) == N || ArgumentError("Number of antenna channels must match")
    p = Pattern(m, x -> abs(s' * x) / norm(x); kwargs...)
    _, idx = findmax(p.values)
    az = p.azs[idx.I[2]]
    el = p.els[idx.I[1]]
    return Spherical(1, az, π/2 - el)
end
