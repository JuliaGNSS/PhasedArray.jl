"""
  plot_pattern(get_steer_vec, reduce_ant_fun = norm, num_az = 360, num_el = 91)

Plots the pattern in polar plot. `reduce_ant_fun` should be the function to reduce the number of antennas to one value.

# Examples
```julia-repl
julia> steer_vec = manifold(0.1904 / 4 * [1 -1 1 -1; 1 1 -1 -1; 0 0 0 0], 1575420e3)
julia> plot_pattern(steer_vec)
```
"""
function plot_pattern(get_steer_vec, reduce_ant_fun = norm, num_az = 360, num_el = 91, fig = figure(), position = (1,1,1))
    azs, els, values = pattern_plotting_data(get_steer_vec, reduce_ant_fun, num_az, num_el)
    ax = draw_polar_axes(fig, position)
    pattern_plot = ax[:pcolormesh](azs, els * 180 / π, values)
    ax[:grid](true, linestyle = "dashed")
    cb = fig[:colorbar](pattern_plot, pad = 0.14)
    cb[:set_label]("Normierte Leistung")
    ax, azs, els, values
end

function pattern_plotting_data(get_steer_vec, reduce_ant_fun = norm, num_az = 360, num_el = 91)
    azs = range(0, stop = 2 * π, length = num_az)
    els = range(0, stop = π / 2, length = num_el)
    values = [reduce_ant_fun(get_steer_vec(Spherical(1.0, az, π / 2 - el))) for el in els, az in azs]
    azs, els, values
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
function plot_pattern_3D(fig, position, get_steer_vec, reduce_ant_fun = norm, num_az = 360, num_el = 181, max_el = num_el - 1)
    X, Y, Z, gains = pattern_3D_plotting_data(get_steer_vec, reduce_ant_fun, num_az, num_el, max_el)
    ax = fig[:add_subplot](position...)
    ax[:plot_surface](X, Y, Z, alpha = 0.4, facecolors = get_cmap("jet")(gains / maximum(gains)), shade = false, linewidth = 0, antialiased = false)
end

function pattern_3D_plotting_data(get_steer_vec, reduce_ant_fun = norm, num_az = 360, num_el = 181, max_el = num_el - 1)
    azs = range(0, stop = 2 * π, length = num_az)
    els = range(0, stop = max_el * π / 180, length = num_el)
    doas_sph = [Spherical(1.0, az, π / 2 - el) for el in els, az in azs]
    doas_cart = [CartesianFromSpherical()(doa_sph) for doa_sph in doas_sph]
    gains = [reduce_ant_fun(get_steer_vec(doa_sph)) for doa_sph in doas_sph]
    scaled_doas_cart = doas_cart .* gains
    X = [doa[1] for doa in scaled_doas_cart]
    Y = [doa[2] for doa in scaled_doas_cart]
    Z = [doa[3] for doa in scaled_doas_cart]
    X, Y, Z, gains
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
function plot_manifold_3D(fig, position, get_steer_vec, ant_pos, reduce_to_real = abs, num_az = 360, num_el = 181, max_el = num_el - 1)
    ax = fig[:add_subplot](position...)
    num_ants = size(ant_pos, 2)
    azs = range(0, stop = 2 * π, length = num_az)
    els = range(0, stop = max_el * π / 180, length = num_el)
    doas_sph = [Spherical(1.0, az, π / 2 - el) for el in els, az in azs]
    doas_cart = [CartesianFromSpherical()(doa_sph) for ant = 1:num_ants, doa_sph in doas_sph]
    gains = [reduce_to_real(get_steer_vec(doa_sph)[ant]) for ant = 1:num_ants, doa_sph in doas_sph]
    scaled_doas_cart = doas_cart .* gains
    X = [doa[1] for doa in scaled_doas_cart]
    Y = [doa[2] for doa in scaled_doas_cart]
    Z = [doa[3] for doa in scaled_doas_cart]
    for ant = 1:num_ants
        ax[:plot_surface](X[ant,:,:] + ant_pos[1,ant], Y[ant,:,:] + ant_pos[2,ant], Z[ant,:,:] + ant_pos[3,ant], alpha = 0.4, facecolors = get_cmap("jet")(gains[ant,:,:] / maximum(gains[ant,:,:])), shade = false, linewidth = 0, antialiased = false)
    end
end

function draw_polar_axes(fig, position)
    ax = fig[:add_subplot](position..., polar="true")
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
