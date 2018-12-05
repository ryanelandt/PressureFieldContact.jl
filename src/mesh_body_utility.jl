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

function Tri_Tet_Intersections.area(mc::MeshCache{Tri,T2}) where {T2}
    area_float = 0.0
    for ind_tri_k = get_ind_tri(mc)
        area_float += area(get_point(mc)[ind_tri_k])
    end
    return area_float
end

function Tri_Tet_Intersections.volume(mc::MeshCache{T1,Tet}) where {T1}
    vol_float = 0.0
    for ind_tet_k = get_ind_tet(mc)
        vol_float += volume(get_point(mc)[ind_tet_k])
    end
    return vol_float
end
