@pyimport matplotlib.animation as anim

"""
# Example
```
using PyPlot
ant_pos = 0.1904 / 4 * [1 -1 1 -1; 1 1 -1 -1; 0 0 0 0]
get_steer_vec = manifold(ant_pos, 1575420e3)
animate_pattern_data = init_animate_pattern_data(get_steer_vec)
animate_pattern(figure(), (1,1,1), animate_pattern_data)
```
"""
function animate_pattern(fig, position, animate_pattern_data)
    azs, els, gains, sat_cn0, sat_doa, jammer_elevation = animate_pattern_data(0.0)
    ax = draw_polar_axes(fig, position)
    pattern_plot = ax[:pcolormesh](azs, els * 180 / π, gains, shading = "gouraud", cmap = get_cmap("jet"))
    ax[:grid](true)
    jammer_plot = ax[:scatter](0.0, 90 - jammer_elevation * 180 / π, c = "r", marker = "x", s = 60)
    sat_plot = ax[:scatter](sat_doa.θ, 90 - sat_doa.ϕ * 180 / π, c = "b", marker = "+", s = 80)
    fig[:colorbar](pattern_plot, label = "Amplification (dB)")
    snr_text = ax[:text](315 * π / 180, 130, @sprintf("CN0: %.1f dB-Hz", sat_cn0))

    function init()
        (pattern_plot, jammer_plot, snr_text)
    end

    # Animate draws the i-th frame, where i starts at i=0 as in Python.
    function animate(i)
        azs, els, gains, sat_cn0, sat_doa, jammer_elevation = animate_pattern_data(i * π / 180)
        pattern_plot[:set_array](vec(gains'))
        jammer_plot[:set_offsets]([i * π / 180, 90 - jammer_elevation * 180 / π])
        snr_text[:set_text](@sprintf("CN0: %.1f dB-Hz", sat_cn0))
        (pattern_plot, jammer_plot, snr_text)
    end

    # Create the animation object by calling the Python function FuncAnimaton
    anim.FuncAnimation(fig, animate, init_func = init, frames = 360, repeat = false, interval = 45)
end

function animate_pattern(fig, position, animate_pattern_data, filename)
    animation = animate_pattern(fig, position, animate_pattern_data)
    animation[:save](filename, bitrate = -1, codec = "libx264", extra_args=["-pix_fmt", "yuv420p"])
end

"""
# Example
```
ant_pos = 0.1904 / 4 * [1 -1 1 -1; 1 1 -1 -1; 0 0 0 0]
get_steer_vec = manifold(ant_pos, 1575420e3)
sat_doa = Spherical(1.0, 10 * π / 180, 75 * π / 180)
init_animate_pattern_data(get_steer_vec, sat_doa, 20, 15 * π / 180)
```
"""
function init_animate_pattern_data(get_steer_vec, sat_doa = Spherical(1.0, 10 * π / 180, 75 * π / 180), jammer_power = 20, jammer_elevation = 15 * π / 180)
    jammer_azimuth -> begin
        jammer_steer_vec = get_steer_vec(Spherical(1.0, jammer_azimuth, jammer_elevation))
        num_ants = length(jammer_steer_vec)
        sat_steer_vec = get_steer_vec(sat_doa)
        noise = complex.(randn(num_ants, 100), randn(num_ants, 100)) / sqrt(2)
        signal = jammer_steer_vec .* 10^(jammer_power / 10) .+ noise
        Rxx = signal * signal' / 100
        prew_filter = Rxx^(-1/2)
        sat_cn0 = 10 * log10(norm(prew_filter * sat_steer_vec)) + 45
        reduce_ant_fun(steer_vec) = 10 * log10(norm(prew_filter * steer_vec))
        azs, els, gains = pattern_plotting_data(get_steer_vec, reduce_ant_fun, 360, 91)
        azs, els, gains, sat_cn0, sat_doa, jammer_elevation
    end
end
