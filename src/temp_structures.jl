
@RigidBodyDynamics.indextype BristleID

struct BristleFriction
    BristleID::BristleID
    tau::Float64
    K_θ::Float64
    K_r::Float64
    function BristleFriction(bristle_ID::BristleID; tau::Float64, K_θ::Float64, K_r::Float64)
        return new(bristle_ID, tau, K_θ, K_r)
    end
end

mutable struct ContactInstructions
    id_tri::MeshID
    id_tet::MeshID
    frac_epsilon::Float64
    frac_linear_weight::Float64
    mu_pair::Float64
    BristleFriction::Union{Nothing,BristleFriction}
    function ContactInstructions(id_tri::MeshID, id_tet::MeshID, frac_epsilon::Float64, frac_linear_weight::Float64, mu_pair::Float64,
        fric_model::Union{Nothing,BristleFriction})

        return new(id_tri, id_tet, frac_epsilon, frac_linear_weight, mu_pair, fric_model)
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
### Volume ###
function add_body_volume_mesh!(ts::TempContactStruct, name::String, point::Vector{SVector{3,Float64}},
    tri_ind::Vector{SVector{3,Int64}}, tet_ind::Vector{SVector{4,Int64}}, strain::Vector{Float64},
    contact_prop::ContactProperties, inertia_prop::InertiaProperties)

    body = add_body_volume!(ts.mechanism, name, point, tet_ind, inertia_prop)
    add_volume_mesh!(ts, body, name, point, tri_ind, tet_ind, strain, contact_prop, inertia_prop)
end
function add_volume_mesh!(ts::TempContactStruct, body::RigidBody{Float64}, name::String, point::Vector{SVector{3,Float64}},
    tri_ind::Vector{SVector{3,Int64}}, tet_ind::Vector{SVector{4,Int64}}, strain::Vector{Float64},
    contact_prop::ContactProperties, inertia_prop::Union{Nothing,InertiaProperties}=nothing)

    mesh = MeshCache(point, name, tri_ind, tet_ind, strain, contact_prop, body, inertia_prop)
    addMesh!(ts, mesh)
end
function add_body_volume!(mechanism::Mechanism, name::String, point::Vector{SVector{3,Float64}}, tet_ind::Vector{SVector{4,Int64}}, inertia_prop::InertiaProperties)
    rho = inertia_prop.rho
    (inertia_prop.thickness == nothing) || error("assumed thickness is something but should be nothing")
    I3, com, mass, mesh_vol = makeInertiaTensor(point, tet_ind, rho)
    return add_body_from_inertia!(mechanism, name, I3, com, mass)
end
### Surface ###
function add_body_surface_mesh!(ts::TempContactStruct, name::String, point::Vector{SVector{3,Float64}},
    tri_ind::Vector{SVector{3,Int64}}, inertia_prop::InertiaProperties)

    body = add_body_surface!(ts.mechanism, name, point, tri_ind, inertia_prop)
    add_surface_mesh!(ts, body, name, point, tri_ind, inertia_prop)
end
function add_surface_mesh!(ts::TempContactStruct, body::RigidBody{Float64}, name::String, point::Vector{SVector{3,Float64}},
        tri_ind::Vector{SVector{3,Int64}}, inertia_prop::Union{Nothing, InertiaProperties}=nothing)

    mesh = MeshCache(point, name, tri_ind, body, inertia_prop)
    addMesh!(ts, mesh)
end
function add_body_surface!(mechanism::Mechanism, name::String, point::Vector{SVector{3,Float64}}, tri_ind::Vector{SVector{3,Int64}}, inertia_prop::InertiaProperties)
    rho = inertia_prop.rho
    thickness = inertia_prop.thickness
    (thickness == nothing) && error("assumed thickness is nothing")
    I3, com, mass, mesh_vol = makeInertiaTensor(point, tri_ind, rho, thickness)
    return add_body_from_inertia!(mechanism, name, I3, com, mass)
end
###
function add_body_from_inertia!(mechanism::Mechanism, name::String, I3, com, mass)
    body_child = newBodyFromInertia(name, I3, com, mass)
    body_parent = root_body(mechanism)
    j_parent_child, x_parent_child = outputJointTransform_ParentChild(body_parent, body_child, SPQuatFloating{Float64}(), SVector{3,Float64}(0,0,0))
    attach!(mechanism, body_parent, body_child, j_parent_child, joint_pose=x_parent_child)
    return body_child
end

function findmesh(ts::MeshCacheDict{MeshCache}, name::String)  # TODO: make this function more elegant
    id = MeshID(-9999)
    for k = keys(ts)
        if ts[k].name == name
            (id == MeshID(-9999)) || error("multiple")
            id = k
        end
    end
    return id
end
