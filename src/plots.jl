"""
  plot_pattern(get_steer_vec, reduce_ant_fun = norm, num_az = 360, num_el = 91)

Plots the pattern in polar plot. `reduce_ant_fun` should be the function to reduce the number of antennas to one value.

# Examples
```julia-repl
julia> steer_vec = manifold(0.1904 / 4 * [1 -1 1 -1; 1 1 -1 -1; 0 0 0 0], 1575420e3)
julia> plot_pattern(steer_vec)
```
"""
function plot_pattern(get_steer_vec, reduce_ant_fun = norm, num_az = 360, num_el = 91)
    azs = linspace(0, 2 * π, num_az)
    els = linspace(0, π / 2, num_el)
    values = [reduce_ant_fun(get_steer_vec(Spherical(1.0, az, π / 2 - el))) for el in els, az in azs]
    figure()
    ax = draw_polar_axes()
    pcolormesh(azs, els * 180 / π, values, cmap = get_cmap("jet"))
    ax[:grid](true)
    colorbar()
end

"""
  plot_pattern_3D(get_steer_vec, reduce_ant_fun = norm, num_az = 360, num_el = 181, max_el = num_el - 1)

Plots the pattern 3D. `reduce_ant_fun` should be the function to reduce the number of antennas to one value.

# Examples
```julia-repl
julia> steer_vec = manifold(0.1904 / 4 * [1 -1 1 -1; 1 1 -1 -1; 0 0 0 0], 1575420e3)
julia> plot_pattern_3D(steer_vec)
```
"""

function plot_pattern_3D(get_steer_vec, reduce_ant_fun = norm, num_az = 360, num_el = 181, max_el = num_el - 1)
    figure()
    azs = linspace(0, 2 * π, num_az)
    els = linspace(0, max_el * π / 180, num_el)
    doas_sph = [Spherical(1.0, az, π / 2 - el) for el in els, az in azs]
    doas_cart = [CartesianFromSpherical()(doa_sph) for doa_sph in doas_sph]
    gains = [reduce_ant_fun(get_steer_vec(doa_sph)) for doa_sph in doas_sph]
    scaled_doas_cart = doas_cart .* gains
    X = [doa[1] for doa in scaled_doas_cart]
    Y = [doa[2] for doa in scaled_doas_cart]
    Z = [doa[3] for doa in scaled_doas_cart]
    plot_surface(X, Y, Z, alpha = 0.4, facecolors = get_cmap("jet")(gains / maximum(gains)), shade = false, linewidth = 0, antialiased = false)
end

"""
  plot_manifold_3D(get_steer_vec, ant_pos, reduce_to_real = abs, num_az = 360, num_el = 181, max_el = num_el - 1)

Plots the manifold for all antenna positions.

# Examples
```julia-repl
julia> steer_vec = manifold(0.1904 / 4 * [1 -1 1 -1; 1 1 -1 -1; 0 0 0 0], 1575420e3)
julia> plot_manifold_3D(steer_vec,[1 -1 1 -1; 1 1 -1 -1; 0 0 0 0])
```
"""

function plot_manifold_3D(get_steer_vec, ant_pos, reduce_to_real = abs, num_az = 360, num_el = 181, max_el = num_el - 1)
    figure()
    num_ants = size(ant_pos, 2)
    azs = linspace(0, 2 * π, num_az)
    els = linspace(0, max_el * π / 180, num_el)
    doas_sph = [Spherical(1.0, az, π / 2 - el) for el in els, az in azs]
    doas_cart = [CartesianFromSpherical()(doa_sph) for ant = 1:num_ants, doa_sph in doas_sph]
    gains = [reduce_to_real(get_steer_vec(doa_sph)[ant]) for ant = 1:num_ants, doa_sph in doas_sph]
    scaled_doas_cart = doas_cart .* gains
    X = [doa[1] for doa in scaled_doas_cart]
    Y = [doa[2] for doa in scaled_doas_cart]
    Z = [doa[3] for doa in scaled_doas_cart]
    for ant = 1:num_ants
        plot_surface(X[ant,:,:] + ant_pos[1,ant], Y[ant,:,:] + ant_pos[2,ant], Z[ant,:,:] + ant_pos[3,ant], alpha = 0.4, facecolors = get_cmap("jet")(gains[ant,:,:] / maximum(gains[ant,:,:])), shade = false, linewidth = 0, antialiased = false)
    end
end

function draw_polar_axes()
    ax = axes(polar="true")
    dtheta=45
    ax[:set_thetagrids](collect(0:dtheta:360-dtheta))
    ax[:set_theta_zero_location]("N")
    ax[:set_theta_direction](-1)

    dr=20
    ax[:set_rgrids](collect(-10:dr:110-dr))
    ax[:set_rlim](0,90)
    ax[:set_rlabel_position](45)
    ax[:set_yticklabels](["","80°","60°","40°","20°",""])
    ax
end
