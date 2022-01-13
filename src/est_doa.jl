"""
$(SIGNATURES)

Calculate the Direction Of Arrival (DOA) with respect to the manifold
`manifold`. The calculated DOA is returned as `Spherical` object, where the
azimuth is counted counter-clockwise from the x axis and the elevation is
counted from the xy-plane towards the zenith.
If the initial azimuth `init_az` and elevation `init_el` are provided, this
will use the hill climbing algorithm.
The hill climbing algorithm should be a lot faster than plain maximum search,
but may only find a local maximum, which might be desired.

# Examples
Using the signal subspace:
```julia-repl
julia> est_doa(manifold, a -> norm(signal_subspace' * a) / norm(a))
```

Using the noise subspace (MUSIC):
```julia-repl
julia> est_doa(manifold, a -> norm(a) / norm(noise_subspace' * a))
```

# Arguments:
- `manifold::AbstractManifold`: Manifold of the antenna array
- `reduction_function`: Function that reduces the manifold vector to a scalar
"""
function est_doa(
    manifold::AbstractManifold,
    reduction_function;
    init_az = nothing,
    init_el = nothing,
    kwargs...
)
    p = Pattern(manifold, reduction_function; kwargs...)
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

This is a convenience function of `est_doa`.
See `est_doa` for more details.
"""
function est_doa_by_noise_subspace(
    manifold::AbstractManifold{N},
    noise_subspace;
    init_az = nothing,
    init_el = nothing,
    kwargs...
) where N
    size(noise_subspace, 1) == N || ArgumentError("Number of antenna channels must match")
    est_doa(
        manifold,
        a -> norm(a) / norm(noise_subspace' * a);
        init_az = init_az,
        init_el = init_el,
        kwargs...
    )
end

"""
$(SIGNATURES)

This is a convenience function of `est_doa`.
See `est_doa` for more details.
"""
function est_doa_by_signal_subspace(
    manifold::AbstractManifold{N},
    signal_subspace;
    init_az = nothing,
    init_el = nothing,
    kwargs...
) where N
    size(signal_subspace, 1) == N || ArgumentError("Number of antenna channels must match")
    est_doa(
        manifold,
        a -> norm(signal_subspace' * a) / norm(a);
        init_az = init_az,
        init_el = init_el,
        kwargs...
    )
end

"""
$(SIGNATURES)

This is a convenience function of `est_doa`.
See `est_doa` for more details.
"""
est_doa_by_music(m, ns; init_az = nothing, init_el = nothing, k...) =
    est_doa_by_noise_subspace(m, ns, init_az = init_az, init_el = init_el, k...)

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