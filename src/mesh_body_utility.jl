function newBodyFromInertia(nameBody::String, I3::SMatrix{3,3,Float64,9}, com::SVector{3,Float64}, m::Float64)
    skew_com = Spatial.vector_to_skew_symmetric(com)
    term = SMatrix{3,3,Float64}(I3 + m * skew_com * transpose(skew_com))
    return RigidBody(SpatialInertia(CartesianFrame3D(nameBody), term, com * m, m))
end

function outputJointTransform_ParentChild(body_parent::RigidBody, body_child::RigidBody, evaluated_joint_type_in, relative_position::SVector{3, Float64})
    j_parent_child = Joint(body_parent.name * "_" * body_child.name, evaluated_joint_type_in)
    x_parent_child = Transform3D(frame_before(j_parent_child), default_frame(body_parent), Quat{Float64}(1,0,0,0), relative_position)
    return j_parent_child, x_parent_child
end
