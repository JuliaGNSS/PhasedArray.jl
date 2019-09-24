abstract type AbstractManifold end

function get_steer_vec(manifold, doa, attitude)
    get_steer_vec(manifold, attitude * doa)
end

function get_steer_vec(manifold, doa::Spherical, attitude)
    get_steer_vec(manifold, sph2cart(doa), attitude)
end
