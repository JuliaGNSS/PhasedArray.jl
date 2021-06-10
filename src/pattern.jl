struct Pattern{
        A<:StepRangeLen,
        E<:StepRangeLen,
        V<:Matrix
    }
    azs::A
    els::E
    values::V
    max_el::Float64
end

struct Pattern3D{T}
    X::Matrix{T}
    Y::Matrix{T}
    Z::Matrix{T}
    gains::Matrix{T}
    max_el::Float64
end

function Pattern(manifold, reduce_ant_fun = norm; num_az = 360, num_el = 91, max_el = π / 2)
    azs = range(0, step = 2π/num_az, length=num_az)
    els = range(0, stop = max_el, length = num_el)
    values = [reduce_ant_fun(get_steer_vec(manifold, Spherical(1.0, az, π / 2 - el))) for el in els, az in azs]
    Pattern(azs, els, values, max_el)
end

@recipe function f(pattern::Pattern;)
    seriestype := :heatmap
    seriescolor --> :viridis
    colorbar_title --> "Magnitude"
    projection --> :polar
    ylims --> (0, pattern.max_el)
    pattern.azs, pattern.els, pattern.values
end

function Pattern3D(manifold, reduce_ant_fun = norm; num_az = 360, num_el = 181, max_el = π)
    azs = range(0, step = 2π/num_az, length=num_az)
    els = range(0, stop = max_el, length = num_el)
    doas_sph = [Spherical(1.0, az, π / 2 - el) for el in els, az in azs]
    doas_cart = map(doa_sph -> CartesianFromSpherical()(doa_sph), doas_sph)
    gains = [reduce_ant_fun(get_steer_vec(manifold, doa_sph)) for doa_sph in doas_sph]
    scaled_doas_cart = doas_cart .* gains
    X = map(doa -> doa[1], scaled_doas_cart)
    Y = map(doa -> doa[2], scaled_doas_cart)
    Z = map(doa -> doa[3], scaled_doas_cart)
    Pattern3D(X, Y, Z, gains, Float64(max_el))
end

@recipe function f(pattern::Pattern3D;)
    seriestype := :surface
    seriescolor --> :viridis
    pattern.X, pattern.Y, pattern.Z
end
