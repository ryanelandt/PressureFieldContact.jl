@RigidBodyDynamics.indextype BristleID

struct BristleFriction
    BristleID::BristleID
    τ::Float64
    K_θ::Float64
    K_r::Float64
    BristleFriction(bristle_ID::BristleID; τ::Float64, K_θ::Float64, K_r::Float64) = new(bristle_ID, τ, K_θ, K_r)
end

mutable struct ContactInstructions
    id_tri::MeshID
    id_tet::MeshID
    frac_ϵ::Float64
    frac_linear_weight::Float64
    μ_pair::Float64
    BristleFriction::Union{Nothing,BristleFriction}
    function ContactInstructions(id_tri::MeshID, id_tet::MeshID, frac_ϵ::Float64, frac_linear_weight::Float64,
        μ_pair::Float64, fric_model::Union{Nothing,BristleFriction})

        return new(id_tri, id_tet, frac_ϵ, frac_linear_weight, μ_pair, fric_model)
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
function add_body_volume_mesh!(ts::TempContactStruct, name::String, h_mesh::HomogenousMesh, tet_mesh::TetMesh,
        inertia_prop::InertiaProperties; joint::JT=SPQuatFloating{Float64}()) where {JT<:JointType}

    body, joint = add_body_volume!(ts.mechanism, name, h_mesh, tet_mesh, inertia_prop, joint=joint)
    add_volume_mesh!(ts, body, name, h_mesh, tet_mesh, inertia_prop)
    return joint
end

function add_volume_mesh!(ts::TempContactStruct, body::RigidBody{Float64}, name::String, h_mesh::HomogenousMesh,
        tet_mesh::TetMesh, inertia_prop::Union{Nothing,InertiaProperties}=nothing)

    mesh = MeshCache(name, h_mesh, tet_mesh, body, inertia_prop)
    addMesh!(ts, mesh)
end

function add_body_volume!(mechanism::Mechanism, name::String, h_mesh::HomogenousMesh, tet_mesh::TetMesh,
        inertia_prop::InertiaProperties; joint::JT=SPQuatFloating{Float64}()) where {JT<:JointType}

    tet_ind = tet_mesh.tet.ind
    point = get_h_mesh_vertices(h_mesh)
    rho = inertia_prop.rho
    (inertia_prop.d == nothing) || error("assumed thickness is something but should be nothing")
    I3, com, mass, mesh_vol = makeInertiaTensor(point, tet_ind, rho)
    return add_body_from_inertia!(mechanism, name, I3, com, mass, joint=joint)
end

### Surface ###
function add_body_surface_mesh!(ts::TempContactStruct, name::String, h_mesh::HomogenousMesh,
        inertia_prop::InertiaProperties; joint::JT=SPQuatFloating{Float64}()) where {JT<:JointType}

    body, joint = add_body_surface!(ts.mechanism, name, h_mesh, inertia_prop, joint=joint)
    add_surface_mesh!(ts, body, name, h_mesh, inertia_prop)
    return joint
end

function add_surface_mesh!(ts::TempContactStruct, body::RigidBody{Float64}, name::String, h_mesh::HomogenousMesh,
        inertia_prop::Union{Nothing, InertiaProperties}=nothing)

    mesh = MeshCache(name, h_mesh, body, inertia_prop)
    addMesh!(ts, mesh)
end

function add_body_surface!(mechanism::Mechanism, name::String, h_mesh::HomogenousMesh, inertia_prop::InertiaProperties;
        joint::JT=SPQuatFloating{Float64}()) where {JT<:JointType}

    point, tri_ind = extract_HomogenousMesh_face_vertices(h_mesh)

    rho = inertia_prop.rho
    d = inertia_prop.d
    (d == nothing) && error("assumed thickness is nothing")
    I3, com, mass, mesh_vol = makeInertiaTensor(point, tri_ind, rho, d)
    return add_body_from_inertia!(mechanism, name, I3, com, mass, joint=joint)
end

function add_body_from_inertia!(mechanism::Mechanism, name::String, I3, com, mass; joint::JT=SPQuatFloating{Float64}()) where {JT<:JointType}
    body_child = newBodyFromInertia(name, I3, com, mass)
    body_parent = root_body(mechanism)
    j_parent_child, x_parent_child = outputJointTransform_ParentChild(body_parent, body_child, joint, SVector{3,Float64}(0,0,0))
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

findMesh(ts::MeshCacheDict{MeshCache}, name::String) = ts[findmesh(ts, name)]
add_pair_rigid_compliant_regularize!(ts::TempContactStruct, name_tri::String, name_tet::String) = add_pair_rigid_compliant!(ts, name_tri, name_tet, nothing)

function add_pair_rigid_compliant!(ts::TempContactStruct, name_tri::String, name_tet::String, friction_model::Union{Nothing,BristleFriction})
    mesh_id_tri = findmesh(ts.MeshCache, name_tri)
    mesh_id_tet = findmesh(ts.MeshCache, name_tet)
    (1 <= mesh_id_tri) || error("invalid tri mesh id $mesh_id_tri")
    (1 <= mesh_id_tet) || error("invalid tet mesh id $mesh_id_tet")
    (mesh_id_tri == mesh_id_tet) && error("tri_mesh and tet_mesh id are the same $mesh_id_tri")
    mesh_cache_tri = ts.MeshCache[mesh_id_tri]
    mesh_cache_tet = ts.MeshCache[mesh_id_tet]
    (mesh_cache_tet.tet == nothing) && error("compliant mesh named $name_tet has no volume mesh")
    mat_tet = mesh_cache_tet.tet.c_prop
    if mesh_cache_tri.tet == nothing
        μ = mat_tet.μ
        frac_ϵ = 1.0
    else
        mat_tri = mesh_cache_tri.c_prop
        μ = calcMutualMu(mat_tri, mat_tet)
        hC_tri = calculateExtrensicCompliance(mat_tri)
        hC_tet = calculateExtrensicCompliance(mat_tet)
        (hC_tet == 0.0) && error("compliance f tet mesh is rigid because its compliance is zero")
        frac_ϵ = hC_tet / (hC_tri + hC_tet)
    end
    (0.0 <= μ <= 3.0) || error("mu our of range")
    frac_linear_weight = 1.0
    new_contact = ContactInstructions(mesh_id_tri, mesh_id_tet, frac_ϵ, frac_linear_weight, μ, friction_model)
    push!(ts.ContactInstructions, new_contact)
    return nothing
end

function add_pair_rigid_compliant_bristle!(ts::TempContactStruct, name_tri::String, name_tet::String; τ::Float64=10.0, K_θ::Union{Nothing,Float64}=nothing, K_r::Union{Nothing, Float64}=nothing)
    (K_θ == nothing) && error("K_θ needs to be given")
    (K_r == nothing) && error("K_r needs to be given")
    bristle_id = BristleID(1 + length(ts.bristle_ids))
    bf = BristleFriction(bristle_id, τ=τ, K_θ=K_θ, K_r=K_r)
    ts.bristle_ids = Base.OneTo(bristle_id)
    return add_pair_rigid_compliant!(ts, name_tri, name_tet, bf)
end

function add_pair_rigid_compliant_bristle_tune_tet!(ts::TempContactStruct, name_tri::String, name_tet::String; τ::Float64=10.0, f_disp::Float64=0.0025, rad_disp::Float64=deg2rad(0.25))
    K_θ, K_r = tune_bristle_stiffness(ts, name_tet, f_disp, rad_disp)
    add_pair_rigid_compliant_bristle!(ts, name_tri, name_tet, τ=τ, K_θ=K_θ, K_r=K_r)
    return nothing
end

function add_pair_rigid_compliant_bristle_tune_tri!(ts::TempContactStruct, name_tri::String, name_tet::String; τ::Float64=10.0, f_disp::Float64=0.0025, rad_disp::Float64=deg2rad(0.25))
    K_θ, K_r = tune_bristle_stiffness(ts, name_tri, f_disp, rad_disp)
    add_pair_rigid_compliant_bristle!(ts, name_tri, name_tet, τ=τ, K_θ=K_θ, K_r=K_r)
    return nothing
end

function tune_bristle_stiffness(ts::TempContactStruct, name::String, f_disp::Float64=0.0025, rad_disp::Float64=deg2rad(0.25))
    K_θ = calc_angular_stiffness(ts, name, rad_disp=rad_disp)
    K_r = calc_linear_stiffness(ts, name, f_disp=f_disp)
    return K_θ, K_r
end

function calc_linear_stiffness(ts::TempContactStruct, name::String; f_disp::Float64=0.0025)
    mesh_id = findmesh(ts.MeshCache, name)
    gravity = norm(ts.mechanism.gravitational_acceleration.v)
    mesh = ts.MeshCache[mesh_id]
    inertia_tensor = bodies(ts.mechanism)[mesh.BodyID].inertia
    F = inertia_tensor.mass * gravity
    char_length = sum(mesh.tri.tree.box.e) / 3  # to avoid importing Statistics
    delta_x = char_length * f_disp
    return F / delta_x
end

function calc_angular_stiffness(ts::TempContactStruct, name::String; rad_disp::Float64=deg2rad(0.25))
    mesh_id = findmesh(ts.MeshCache, name)
    mesh = ts.MeshCache[mesh_id]
    inertia_tensor = bodies(ts.mechanism)[mesh.BodyID].inertia
    avg_inertia = sum(svd(inertia_tensor.moment).S) / 3  # to avoid importing Statistics
    return avg_inertia / rad_disp
end
