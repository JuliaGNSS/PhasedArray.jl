# The golden spiral method
# http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
# min_phi defines the minimum elevation. 90° is the zenith, 0° horizon and -90° the nadir
function create_points_on_sphere(num_points, min_phi = -20 * π / 180, ϵ = 0.36)
    map(0:num_points - 1) do index
        phi = acos(1 - (1 - sin(min_phi)) * (index + ϵ) / (num_points - 1 + 2 * ϵ))
        theta = π * (1 + sqrt(5)) * index
        SVector(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi))
    end
end