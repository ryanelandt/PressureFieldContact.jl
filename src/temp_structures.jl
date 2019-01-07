@RigidBodyDynamics.indextype BristleID

struct Bristle
    BristleID::BristleID
    τ::Float64
    k̄::Float64
    fric_pro::Float64
    Bristle(bristle_ID::BristleID; τ::Float64, k̄::Float64, fric_pro::Float64=2.0) = new(bristle_ID, τ, k̄, fric_pro)
end

struct Regularized
    v_tol⁻¹::Float64
    Regularized(v_tol) = new(1 / v_tol)
end

mutable struct ContactInstructions
    id_1::MeshID
    id_2::MeshID
    mutual_compliance::Bool
    FrictionModel::Union{Regularized,Bristle}
    μ::Float64
    χ::Float64
    function ContactInstructions(id_tri::MeshID, id_tet::MeshID, mutual_compliance::Bool,
            fric_model::Union{Regularized,Bristle}; μ::Float64, χ::Float64)

        # (0.001 <= χ <= 5.0) || error("hc_velocity_damping in unexpected range.")
        (0.0 <= μ <= 3.0) || error("mu in unexpected range.")

        return new(id_tri, id_tet, mutual_compliance, fric_model, μ, χ)
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

function add_pair_rigid_compliant_regularize!(ts::TempContactStruct, name_tri::String, name_tet::String;
        μ::Union{Nothing,Float64}=nothing, χ::Union{Nothing,Float64}=nothing, v_tol::Union{Nothing,Float64}=nothing)
    #
    if v_tol == nothing
        @warn("unspecified v_tol replaced with 0.25")
        v_tol = 0.25
    end
    regularized = Regularized(v_tol)
    return add_pair_rigid_compliant!(ts, name_tri, name_tet, regularized, μ=μ, χ=χ)
end

function add_pair_rigid_compliant!(ts::TempContactStruct, name_1::String, name_2::String,
        friction_model::Union{Regularized,Bristle}; μ::Union{Nothing,Float64}=nothing,
        χ::Union{Nothing,Float64}=nothing)

    mesh_id_1 = findmesh(ts.MeshCache, name_1)
    mesh_id_2 = findmesh(ts.MeshCache, name_2)
    (1 <= mesh_id_1) || error("invalid 1 mesh id $mesh_id_1")
    (1 <= mesh_id_2) || error("invalid 2 mesh id $mesh_id_2")
    (mesh_id_1 == mesh_id_2) && error("1_mesh and tet_mesh id are the same $mesh_id_1")
    mesh_cache_1 = ts.MeshCache[mesh_id_1]
    mesh_cache_2 = ts.MeshCache[mesh_id_2]
    is_compliant_1 = is_compliant(mesh_cache_1)
    is_compliant_2 = is_compliant(mesh_cache_2)
    is_compliant_1 || is_compliant_2 || error("neither mesh is compliant")
    if is_compliant_1
        mesh_id_1, mesh_id_2 = mesh_id_2, mesh_id_1
    end
    if μ == nothing
        @warn("unspecified μ replaced with 0.3")
        μ = 0.3
    end
    if χ == nothing
        @warn("unspecified χ replaced with 0.5")
        χ = 0.5
    end
    mutual_compliance = is_compliant_1 && is_compliant_2
    push!(ts.ContactInstructions, ContactInstructions(mesh_id_1, mesh_id_2, mutual_compliance, friction_model, μ=μ, χ=χ))
    return nothing
end

function add_pair_rigid_compliant_bristle!(ts::TempContactStruct, name_tri::String, name_tet::String; τ::Float64=30.0,
        k̄=1.0e4, fric_pro=2.0, μ::Union{Nothing,Float64}=nothing, χ::Union{Nothing,Float64}=nothing)

    bristle_id = BristleID(1 + length(ts.bristle_ids))
    bf = Bristle(bristle_id, τ=τ, k̄=k̄, fric_pro=fric_pro)
    ts.bristle_ids = Base.OneTo(bristle_id)
    return add_pair_rigid_compliant!(ts, name_tri, name_tet, bf, μ=μ, χ=χ)
end
