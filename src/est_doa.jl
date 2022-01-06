"""
$(SIGNATURES)

Calculate the direction of arrival of the signal `s` with respect to the 
manifold `m`. The calculated DOA is returned as `Spherical` object, where the
azimuth is counted counter-clockwise from the x axis and the elevation is
counted from the xy-plane towards the zenith.
If the initial azimuth `init_az` and elevation `init_el` is provided, this
will use the hill climbing algorithm.
The hill climbing algorithm should be a lot faster than plain maximum search,
but may only find a local maximum, which might be desired.

# Arguments:
- `m::AbstractManifold{N}`: Manifold of the antenna array
- `s::AbstractVector`: Signal subspace for which the doa is estimated
"""
function est_doa_by_music(
    m::AbstractManifold{N},
    s::AbstractVector;
    init_az = nothing,
    init_el = nothing,
    kwargs...
) where N
    length(s) == N || ArgumentError("Number of antenna channels must match")
    p = Pattern(m, x -> abs(s' * x) / norm(x); kwargs...)
    if isnothing(init_az) || isnothing(init_el) 
        _, idx = findmax(p.values)
    else
        init_az_idx = floor(Int, mod2pi(init_az) / Float64(p.azs.step) + 1)
        init_el_idx = floor(Int, (π / 2 - init_el) / Float64(p.els.step) + 1)
        _, idx = hill_climbing(p.values, init_el_idx, init_az_idx)
    end
    az = p.azs[idx.I[2]]
    el = p.els[idx.I[1]]
    return Spherical(1, az, π / 2 - el)
end

"""
$(SIGNATURES)

Multi-dimensionial hill climbing. 

# Arguments:
- `values::AbstractArray`: Some multi-dimensional array
- `inits`: Starting indices
"""
function hill_climbing(values, inits...)
    ndims(values) == length(inits) || ArgumentError("Number of inits does not fit dimensions of values")
    indices = CartesianIndices(values)
    lower_index_bounds, upper_index_bounds = first(indices), last(indices)
    step_size = oneunit(lower_index_bounds)
    current_index = CartesianIndex(inits...)
    current_value = values[current_index]
    while true
        next_index = current_index
        for J in max(lower_index_bounds, current_index - step_size):min(upper_index_bounds, current_index + step_size)
            if values[J] > current_value
                next_index = J
                current_value = values[J]
            end
        end
        if current_index == next_index
            return current_value, current_index
        else
            current_index = next_index
        end
    end
end