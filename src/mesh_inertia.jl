makeInertiaTensor(point::Vector{SVector{3,Float64}}, vec_vol_ind::Vector{SVector{3, Int64}}, rho::Float64, thickness::Float64) = makeInertiaTensor(getTriQuadRule(3), point, vec_vol_ind, rho, thickness)
makeInertiaTensor(point::Vector{SVector{3,Float64}}, vec_vol_ind::Vector{SVector{4, Int64}}, rho::Float64, ::Float64) = makeInertiaTensor(getTetQuadRule(4), point, vec_vol_ind, rho)

function makeInertiaTensor(quad_rule::TriTetQuadRule{N_zeta, NQ}, point::Vector{SVector{3,Float64}}, vec_vol_ind::Vector{SVector{N_zeta, Int64}}, rho::Float64, thickness=NaN) where {N_zeta, NQ}
    com, mesh_vol = centroidVolumeCombo(point, vec_vol_ind, thickness)
    tensor_I = zeros(3, 3)
    eye3 = SMatrix{3,3}(1.0I)
    for ind_k = vec_vol_ind
        points_simplex_k = point[ind_k]
        A = asMat(points_simplex_k)
        v_tet = equiv_volume(points_simplex_k, thickness)
        for k_quad_point = 1:NQ
            zeta_quad = quad_rule.zeta[k_quad_point]
            r = A * zeta_quad - com
            raw_tensor = eye3 * dot(r, r) - (r * transpose(r))
            the_mass = rho * quad_rule.w[k_quad_point] * v_tet
            tensor_I += the_mass * raw_tensor
        end
    end
    return tensor_I, com, mesh_vol * rho, mesh_vol
end

function centroidVolumeCombo(point::Vector{SVector{3,Float64}}, vec_vol_ind::Vector{SVector{N, Int64}}, thickness=NaN) where {N}
    v_cum, c_cum = 0.0, 0.0
    for ind_k = vec_vol_ind
        points_simplex_k = point[ind_k]
        v = equiv_volume(points_simplex_k, thickness)
        v_cum += v
        c_cum += v * centroid(points_simplex_k)
    end
    return c_cum / v_cum, v_cum
end

equiv_volume(v::SVector{4,SVector{3,Float64}},::Float64) = volume(v)
equiv_volume(v::SVector{3,SVector{3,Float64}}, h::Float64) = area(v) * h
