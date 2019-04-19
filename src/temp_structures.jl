@RigidBodyDynamics.indextype BristleID

struct Bristle
    BristleID::BristleID
    τ::Float64
    k̄::Float64
    Bristle(bristle_ID::BristleID; τ::Float64, k̄::Float64) = new(bristle_ID, τ, k̄)
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

        (0.0 <= μ <= 3.0) || error("mu in unexpected range.")
        return new(id_tri, id_tet, mutual_compliance, fric_model, μ, χ)
    end
end

mutable struct TempContactStruct
    is_aabb::Bool
    mechanism::Mechanism
    mesh_ids::Base.OneTo{MeshID}
    bristle_ids::Base.OneTo{BristleID}
    MeshCache::RigidBodyDynamics.CustomCollections.CacheIndexDict{MeshID,Base.OneTo{MeshID},MeshCache}
    ContactInstructions::Vector{ContactInstructions}
    function TempContactStruct(mechanism::Mechanism, is_aabb::Bool=false)
        bristle_ids = Base.OneTo(BristleID(0))
        mesh_ids = Base.OneTo(MeshID(0))
        mesh_cache = MeshCacheDict{MeshCache}(mesh_ids)
        vec_ins = Vector{ContactInstructions}()
        return new(is_aabb, mechanism, mesh_ids, bristle_ids, mesh_cache, vec_ins)
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

function add_body_contact!(ts::TempContactStruct, name::String, e_mesh::eMesh;
        i_prop::InertiaProperties,
        c_prop::Union{Nothing,ContactProperties}=nothing,
        body::Union{RigidBody{Float64},Nothing}=nothing,
        joint::JT=SPQuatFloating{Float64}(),
        dh::basic_dh=one(basic_dh{Float64})) where {JT<:JointType}

    nt = add_body!(ts, name, e_mesh, i_prop=i_prop, body=body, joint=joint, dh=dh)
    nt_new_contact = add_contact!(ts, name, e_mesh, c_prop=c_prop, body=nt.body)
    return NamedTuple{(:body, :joint, :id)}((nt.body, nt.joint, nt_new_contact.id))
end

function make_eTree_obb(eM_box::eMesh{T1,T2}, c_prop::Union{Nothing,ContactProperties}) where {T1,T2}

    xor(c_prop == nothing, T2 == Nothing) && error("Attempting to use nothing as the ContartProperties for a Tet mesh")

    e_tree = eTree(eM_box, c_prop)

    if T1 != Nothing
        all_obb_tri = [fit_tri_obb(eM_box, k) for k = 1:n_tri(eM_box)]
        obb_tri = obb_tree_from_aabb(e_tree.tri, all_obb_tri)
    else
        obb_tri = nothing
    end
    if T2 != Nothing
        all_obb_tet = [fit_tet_obb(eM_box, k) for k = 1:n_tet(eM_box)]
        obb_tet = obb_tree_from_aabb(e_tree.tet, all_obb_tet)
    else
        obb_tet = nothing
    end

    return eTree(obb_tri, obb_tet, c_prop)
end

function add_contact!(ts::TempContactStruct, name::String, e_mesh::eMesh;
        c_prop=c_prop::Union{Nothing,ContactProperties}=nothing,
        body::Union{RigidBody{Float64},Nothing}=nothing)

    body = return_body_never_nothing(ts.mechanism, body)
    if ts.is_aabb
        e_tree = eTree(e_mesh, c_prop)
    else
        e_tree = make_eTree_obb(e_mesh, c_prop)
    end
    mesh = MeshCache(name, e_mesh, e_tree, body)
    addMesh!(ts, mesh)
    return NamedTuple{(:id,)}((find_mesh_id(ts, mesh),))
end

function add_body!(ts::TempContactStruct, name::String, e_mesh::eMesh; i_prop::InertiaProperties,
        body::Union{RigidBody{Float64},Nothing}=nothing, joint::JT=SPQuatFloating{Float64}(),
        dh::basic_dh=one(basic_dh{Float64})) where {JT<:JointType}

    mesh_inertia_info = makeInertiaInfo(e_mesh, i_prop)
    return add_body_from_inertia!(ts.mechanism, name, mesh_inertia_info, joint=joint, body=body, dh=dh)
end

return_body_never_nothing(mechanism::Mechanism, body::Nothing) = root_body(mechanism)
return_body_never_nothing(mechanism::Mechanism, body::RigidBody{Float64}) = body

function add_body_from_inertia!(mechanism::Mechanism, name::String, mesh_inertia_info::MeshInertiaInfo;
        joint::JT=SPQuatFloating{Float64}(), body::Union{RigidBody{Float64},Nothing}=nothing,
        dh::basic_dh{Float64}=one(basic_dh{Float64})) where {JT<:JointType}

    body_parent = return_body_never_nothing(mechanism, body)
    body_child = newBodyFromInertia(name, mesh_inertia_info)
    j_parent_child, x_parent_child = outputJointTransform_ParentChild(body_parent, body_child, joint, dh)
    attach!(mechanism, body_parent, body_child, j_parent_child, joint_pose=x_parent_child)
    return NamedTuple{(:body, :joint)}((body_child, j_parent_child))
end

default_χ() = 0.5
default_μ() = 0.3

function add_friction_regularize!(ts::TempContactStruct, mesh_id_1::MeshID, mesh_id_2::MeshID;
        μ::Float64=default_μ(), χ::Float64=default_χ(), v_tol::Float64=0.01)

    regularized = Regularized(v_tol)
    return add_friction!(ts, mesh_id_1, mesh_id_2, regularized, μ=μ, χ=χ)
end

function add_friction!(ts::TempContactStruct, mesh_id_1::MeshID, mesh_id_c::MeshID,
        friction_model::Union{Regularized,Bristle}; μ::Float64, χ::Float64)

    mesh_1 = ts.MeshCache[mesh_id_1]
    mesh_c = ts.MeshCache[mesh_id_c]
    (mesh_1 == mesh_c) && error("mesh_1 and mesh_c are the same")
    is_compliant_1 = is_compliant(mesh_1)
    is_compliant_c = is_compliant(mesh_c)
    is_compliant_1 || is_compliant_c || error("neither mesh is compliant")
    if is_compliant_1
        mesh_id_1, mesh_id_c = mesh_id_c, mesh_id_1
    end
    mutual_compliance = is_compliant_1 && is_compliant_c
    push!(ts.ContactInstructions, ContactInstructions(mesh_id_1, mesh_id_c, mutual_compliance, friction_model, μ=μ, χ=χ))
    return nothing
end

function add_friction_bristle!(ts::TempContactStruct, mesh_id_1::MeshID, mesh_id_c::MeshID;
        τ::Float64=0.05, k̄=1.0e4, μ::Float64=default_μ(), χ::Float64=default_χ())

    isa(μ, Nothing) || (0 < μ) || error("μ cannot be 0 for bristle friction")
    bristle_id = BristleID(1 + length(ts.bristle_ids))
    bf = Bristle(bristle_id, τ=τ, k̄=k̄)
    ts.bristle_ids = Base.OneTo(bristle_id)
    return add_friction!(ts, mesh_id_1, mesh_id_c, bf, μ=μ, χ=χ)
end
