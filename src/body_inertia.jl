
function newBodyFromInertia(nameBody::String, mesh_inertia_info::MeshInertiaInfo)
    com = mesh_inertia_info.com
    m = mesh_inertia_info.mass
    I3 = mesh_inertia_info.tensor_I
    skew_com = Spatial.vector_to_skew_symmetric(com)
    term = SMatrix{3,3,Float64}(I3 + m * skew_com * transpose(skew_com))
    return RigidBody(SpatialInertia(CartesianFrame3D(nameBody), term, com * m, m))
end

function outputJointTransform_ParentChild(body_parent::RigidBody, body_child::RigidBody, evaluated_joint_type_in,
        dh::basic_dh{Float64}=one(basic_dh{Float64}) )

    rot, trans = dh_R_t(dh)
    rot = RotMatrix{3,Float64}(rot)
    j_parent_child = Joint(body_parent.name * "_" * body_child.name, evaluated_joint_type_in)
    x_parent_child = Transform3D(frame_before(j_parent_child), default_frame(body_parent), rot, trans)
    return j_parent_child, x_parent_child
end


function makeInertiaInfo(e_mesh::eMesh{Tri,Nothing}, i_prop::InertiaProperties{Tri})
    return makeInertiaTensor(getTriQuadRule(3), e_mesh.point, e_mesh.tri, i_prop)
end

function makeInertiaInfo(e_mesh::eMesh{Nothing,Tet}, i_prop::InertiaProperties{Tet})
    return makeInertiaTensor(getTetQuadRule(4), e_mesh.point, e_mesh.tet, i_prop)
end

function makeInertiaTensor(quad_rule::TriTetQuadRule{N_ζ,NQ}, point::Vector{SVector{3,Float64}},
        vec_vol_ind::Vector{SVector{N_ζ,Int64}}, i_prop::InertiaProperties) where {N_ζ, NQ}

    d = i_prop.d
    rho = i_prop.rho

    com, mesh_vol = centroidVolumeCombo(point, vec_vol_ind, d)
    tensor_I = zeros(3, 3)
    eye3 = SMatrix{3,3}(1.0I)
    for ind_k = vec_vol_ind
        points_simplex_k = point[ind_k]
        A = asMat(points_simplex_k)
        v_tet = equiv_volume(points_simplex_k, d)
        for k_quad_point = 1:NQ
            ζ_quad = quad_rule.zeta[k_quad_point]
            r = A * ζ_quad - com
            raw_tensor = eye3 * dot(r, r) - (r * transpose(r))
            the_mass = rho * quad_rule.w[k_quad_point] * v_tet
            tensor_I += the_mass * raw_tensor
        end
    end
    return MeshInertiaInfo(tensor_I, com, mesh_vol * rho, mesh_vol)
end

function centroidVolumeCombo(point::Vector{SVector{3,Float64}}, vec_vol_ind::Vector{SVector{N,Int64}},
        d::Union{Nothing,Float64}) where {N}

    v_cum, c_cum = 0.0, 0.0
    for ind_k = vec_vol_ind
        points_simplex_k = point[ind_k]
        v = equiv_volume(points_simplex_k, d)
        v_cum += v
        c_cum += v * centroid(points_simplex_k)
    end
    return c_cum / v_cum, v_cum
end

equiv_volume(v::SVector{4,SVector{3,Float64}}, ::Nothing) = volume(v)
equiv_volume(v::SVector{3,SVector{3,Float64}}, h::Float64) = area(v) * h
