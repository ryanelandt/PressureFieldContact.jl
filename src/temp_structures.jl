@RigidBodyDynamics.indextype BristleID

struct BristleFriction
    BristleID::BristleID
    τ::Float64
    k̄::Float64
    BristleFriction(bristle_ID::BristleID; τ::Float64, k̄::Float64) = new(bristle_ID, τ, k̄)
end

mutable struct ContactInstructions
    id_1::MeshID  # tri or tet
    id_2::MeshID
    mutual_compliance::Bool
    μ_pair::Float64
    BristleFriction::Union{Nothing,BristleFriction}
    function ContactInstructions(id_tri::MeshID, id_tet::MeshID, mutual_compliance::Bool, μ_pair::Float64,
            fric_model::Union{Nothing,BristleFriction})

        return new(id_tri, id_tet, mutual_compliance, μ_pair, fric_model)
    end
end

mutable struct TempContactStruct
    mechanism::Mechanism
    mesh_ids::Base.OneTo{MeshID}
    bristle_ids::Base.OneTo{BristleID}
    MeshCache::RigidBodyDynamics.CustomCollections.CacheIndexDict{MeshID,Base.OneTo{MeshID},MeshCache}
    ContactInstructions::Vector{ContactInstructions}
    function TempContactStruct(mechanism::Mechanism)
        bristle_ids = Base.OneTo(BristleID(0))
        mesh_ids = Base.OneTo(MeshID(0))
        mesh_cache = MeshCacheDict{MeshCache}(mesh_ids)
        vec_ins = Vector{ContactInstructions}()
        return new(mechanism, mesh_ids, bristle_ids, mesh_cache, vec_ins)
    end
end

function addMesh!(ts::TempContactStruct, mesh::MeshCache)
    mesh_ids_old = ts.mesh_ids
    mesh_ids_new = Base.OneTo(MeshID(length(ts.mesh_ids) + 1))
    mesh_cache = MeshCacheDict{MeshCache}(mesh_ids_new)
    for id = mesh_ids_old
        mesh_cache[id] = ts.MeshCache[id]
    end
    mesh_cache[mesh_ids_new[end]] = mesh
    ts.MeshCache = mesh_cache
    ts.mesh_ids = mesh_ids_new
    return nothing
end

function add_body_contact!(ts::TempContactStruct, name::String, e_mesh::eMesh,
        c_prop::Union{Nothing,ContactProperties}, i_prop::InertiaProperties;
        body_parent::Union{RigidBody{Float64},Nothing}=nothing,
        joint_type::JT=SPQuatFloating{Float64}(), dh::basic_dh=one(basic_dh{Float64})) where {JT<:JointType}
    #
    body, joint = add_body!(ts, name, e_mesh, i_prop, body_parent=body_parent, joint_type=joint_type, dh=dh)
    add_contact!(ts, name, e_mesh, c_prop, body=body, dh=dh)
    return body, joint
end

function add_contact!(ts::TempContactStruct, name::String, e_mesh::eMesh,
        c_prop::Union{Nothing,ContactProperties}; body::Union{RigidBody{Float64},Nothing}=nothing,
        dh::basic_dh=one(basic_dh{Float64}))
    #
    body = return_body_never_nothing(ts.mechanism, body)
    e_tree = eTree(e_mesh, c_prop)
    mesh = MeshCache(name, e_mesh, e_tree, body)
    addMesh!(ts, mesh)
    return nothing
end

function add_body!(ts::TempContactStruct, name::String, e_mesh::eMesh, i_prop::InertiaProperties;
        body_parent::Union{RigidBody{Float64},Nothing}=nothing, joint_type::JT=SPQuatFloating{Float64}(),
        dh::basic_dh=one(basic_dh{Float64})) where {JT<:JointType}
    #
    mesh_inertia_info = makeInertiaInfo(e_mesh, i_prop)
    return add_body_from_inertia!(ts.mechanism, name, mesh_inertia_info, joint=joint_type, body_parent=body_parent, dh=dh)
end

return_body_never_nothing(mechanism::Mechanism, body::Nothing) = root_body(mechanism)
return_body_never_nothing(mechanism::Mechanism, body::RigidBody{Float64}) = body

function add_body_from_inertia!(mechanism::Mechanism, name::String, mesh_inertia_info::MeshInertiaInfo;
        joint::JT=SPQuatFloating{Float64}(), body_parent::Union{RigidBody{Float64},Nothing}=nothing,
        dh::basic_dh{Float64}=one(basic_dh{Float64})) where {JT<:JointType}

    body_parent = return_body_never_nothing(mechanism, body_parent)
    body_child = newBodyFromInertia(name, mesh_inertia_info)
    j_parent_child, x_parent_child = outputJointTransform_ParentChild(body_parent, body_child, joint, dh)
    attach!(mechanism, body_parent, body_child, j_parent_child, joint_pose=x_parent_child)
    return body_child, j_parent_child
end

function findmesh(ts::MeshCacheDict{MeshCache}, name::String)  # TODO: make this function more elegant
    id = MeshID(-9999)
    for k = keys(ts)
        if ts[k].name == name
            (id == MeshID(-9999)) || error("multiple")
            id = k
        end
    end
    (id == MeshID(-9999)) && error("no mesh found by name: $name")
    return id
end

function add_pair_rigid_compliant_regularize!(ts::TempContactStruct, name_tri::String, name_tet::String)
    return add_pair_rigid_compliant!(ts, name_tri, name_tet, nothing)
end

function add_pair_rigid_compliant!(ts::TempContactStruct, name_1::String, name_2::String,
        friction_model::Union{Nothing,BristleFriction})

    mesh_id_1 = findmesh(ts.MeshCache, name_1)
    mesh_id_2 = findmesh(ts.MeshCache, name_2)
    (1 <= mesh_id_1) || error("invalid 1 mesh id $mesh_id_1")
    (1 <= mesh_id_2) || error("invalid 2 mesh id $mesh_id_2")
    (mesh_id_1 == mesh_id_2) && error("1_mesh and tet_mesh id are the same $mesh_id_1")
    mesh_cache_1 = ts.MeshCache[mesh_id_1]
    mesh_cache_2 = ts.MeshCache[mesh_id_2]
    is_compliant_1 = is_compliant(mesh_cache_1)
    is_compliant_2 = is_compliant(mesh_cache_2)
    !is_compliant_1 && !is_compliant_2 && error("neither mesh is compliant")
    if is_compliant_1 && !is_compliant_2
        return add_pair_rigid_compliant!(ts, name_2, name_1, friction_model)
    else
        μ = calc_mutual_μ(mesh_cache_1, mesh_cache_2)
        mutual_compliance = is_compliant_1 && is_compliant_2
        new_contact = ContactInstructions(mesh_id_1, mesh_id_2, mutual_compliance, μ, friction_model)
        push!(ts.ContactInstructions, new_contact)
        return nothing
    end
end

function add_pair_rigid_compliant_bristle!(ts::TempContactStruct, name_tri::String, name_tet::String; τ::Float64=30.0,
        k̄=1.0e4)

    bristle_id = BristleID(1 + length(ts.bristle_ids))
    bf = BristleFriction(bristle_id, τ=τ, k̄=k̄)  # , K_θ=K_θ, K_r=K_r)
    ts.bristle_ids = Base.OneTo(bristle_id)
    return add_pair_rigid_compliant!(ts, name_tri, name_tet, bf)
end
