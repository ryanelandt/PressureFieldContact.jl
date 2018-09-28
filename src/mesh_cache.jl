struct SimplexTree{N}
    tree::bin_BB_Tree{AABB}
    ind::Vector{SVector{N,Int64}}
    function SimplexTree(point::Vector{SVector{3,Float64}}, ind::Vector{SVector{3,Int64}})
        tree = triTetMeshToTreeAABB(point, ind)
        return new{3}(tree, ind)
    end
    function SimplexTree(point::Vector{SVector{3,Float64}}, ind::Vector{SVector{4,Int64}})
        tree = triTetMeshToTreeAABB(point, ind)
        return new{4}(tree, ind)
    end
end

struct TetMesh
    tet::SimplexTree{4}
    strain::Vector{Float64}
    contact_prop::ContactProperties
    function TetMesh(point::Vector{SVector{3,Float64}}, tet_ind::Vector{SVector{4,Int64}}, strain::Vector{Float64}, contact_prop::ContactProperties)
        tet_simp_tree = SimplexTree(point, tet_ind)
        return new(tet_simp_tree, strain, contact_prop)
    end
end

struct MeshCache
    point::Vector{SVector{3,Float64}}
    name::String
    BodyID::Union{Nothing,BodyID}
    FrameID::CartesianFrame3D
    InertiaProperties::Union{Nothing, InertiaProperties}
    tri::SimplexTree{3}
    tet::Union{Nothing, TetMesh}

    function MeshCache(point::Vector{SVector{3,Float64}}, name::String, tri_ind::Vector{SVector{3,Int64}}, tet_ind::Vector{SVector{4,Int64}},
        strain::Vector{Float64}, contact_prop::ContactProperties, body::RigidBody{Float64}, inertia_prop::Union{Nothing, InertiaProperties}=nothing)

        tri_simp_tree = SimplexTree(point, tri_ind)
        tet_mesh = TetMesh(point, tet_ind, strain, contact_prop)
        return new(point, name, BodyID(body), default_frame(body), inertia_prop, tri_simp_tree, tet_mesh)
    end

    function MeshCache(point::Vector{SVector{3,Float64}}, name::String, tri_ind::Vector{SVector{3,Int64}}, body::RigidBody{Float64}, inertia_prop::Union{Nothing, InertiaProperties}=nothing)
        tri_simp_tree = SimplexTree(point, tri_ind)
        return new(point, name, BodyID(body), default_frame(body), inertia_prop, tri_simp_tree, nothing)
    end
end

function addBodyMeshCache(point::Vector{SVector{3,Float64}}, name::String, tri_ind::Vector{SVector{3,Int64}}, tet_ind::Vector{SVector{4,Int64}},
    strain::Vector{Float64}, contact_prop::ContactProperties, inertia_prop::InertiaProperties, mechanism::Mechanism)

    body = addBodyVolumeMesh!(mechanism, name, point, tet_ind, inertia_prop)
    return MeshCache(point, name, tri_ind, tet_ind, strain, contact_prop, body, inertia_prop)
end

function addBodyMeshCache(point::Vector{SVector{3,Float64}}, name::String, tri_ind::Vector{SVector{3,Int64}}, inertia_prop::InertiaProperties, mechanism::Mechanism)
    body = addBodySurfaceMesh!(mechanism, name, point, tri_ind, inertia_prop)
    return MeshCache(point, name, tri_ind, body, inertia_prop)
end

function addBodyVolumeMesh!(mechanism::Mechanism, name::String, point::Vector{SVector{3,Float64}}, tet_ind::Vector{SVector{4,Int64}}, inertia_prop::InertiaProperties)
    rho = inertia_prop.rho
    (inertia_prop.thickness == nothing) || error("assumed thickness is something but should be nothing")
    I3, com, mass, mesh_vol = makeInertiaTensor(point, tet_ind, rho)
    return addBodyMesh!(mechanism, name, I3, com, mass)
end

function addBodySurfaceMesh!(mechanism::Mechanism, name::String, point::Vector{SVector{3,Float64}}, tri_ind::Vector{SVector{3,Int64}}, inertia_prop::InertiaProperties)
    rho = inertia_prop.rho
    thickness = inertia_prop.thickness
    (thickness == nothing) && error("assumed thickness is nothing")
    I3, com, mass, mesh_vol = makeInertiaTensor(point, tri_ind, rho, thickness)
    return addBodyMesh!(mechanism, name, I3, com, mass)
end

function addBodyMesh!(mechanism::Mechanism, name::String, I3, com, mass)
    body_child = newBodyFromInertia(name, I3, com, mass)
    body_parent = root_body(mechanism)
    j_parent_child, x_parent_child = outputJointTransform_ParentChild(body_parent, body_child, SPQuatFloating{Float64}(), SVector{3,Float64}(0,0,0))
    attach!(mechanism, body_parent, body_child, j_parent_child, joint_pose=x_parent_child)
    return body_child
end

function asHomogenousMesh(meshCache::MeshCache)
    vec_Face = Face{3, Int32}.(meshCache.tri_simp_tree.ind)  # TODO: consider changing this to Int32
    vec_Point = Point{3, Float32}.(meshCache.point)  # TODO: consider changing this to Float32
    return HomogenousMesh(vec_Point, vec_Face)
end

@RigidBodyDynamics.indextype MeshID
const MeshDict{V} = RigidBodyDynamics.IndexDict{MeshID, Base.OneTo{MeshID}, V}
const MeshCacheDict{V} = RigidBodyDynamics.CacheIndexDict{MeshID, Base.OneTo{MeshID}, V}
Base.@propagate_inbounds Base.getindex(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, key::RigidBody) = d[BodyID(key)]
Base.@propagate_inbounds Base.setindex!(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, value, key::RigidBody) = d[BodyID(key)] = value

# new_mesh_cache() = MeshCacheDict{MeshCache}(Base.OneTo(MeshID(0)))
# function remake_mesh_cache(mesh_cache_old::MeshCacheDict{MeshCache}, name::String, raw::RawMeshCache, body::RigidBody{Float64})
#     mesh_n_old = length(mesh_cache_old)
#     mesh_n_new = mesh_n_old + 1
#     mesh_ids_old = Base.OneTo(MeshID(mesh_n_old))
#     mesh_ids_new = Base.OneTo(MeshID(mesh_n_new))
#     mesh_cache_new = MeshCacheDict{MeshCache}(mesh_ids_new)
#     for id = mesh_ids_old
#         mesh_cache_new[id] = mesh_cache_old[id]
#     end
#     mesh_cache_new[MeshID(mesh_n_new)] = MeshCache(name, raw, body)
#     return mesh_cache_new
# end

# struct RawMeshCache
#     point::Vector{SVector{3,Float64}}
#     tri_ind::Vector{SVector{3,Int64}}
#     tri_tree::bin_BB_Tree{AABB}
#     tet_ind::Union{Nothing,Vector{SVector{4,Int64}}}
#     tet_tree::Union{Nothing,bin_BB_Tree{AABB}}
#     strain::Union{Nothing,Vector{Float64}}
#     material::MaterialProperties
#     function RawMeshCache(point::Vector{SVector{3,Float64}}, tri_ind::Vector{SVector{3,Int64}}, material::MaterialProperties=MaterialProperties())
#         tri_tree = triTetMeshToTreeAABB(point, tri_ind)
#         return new(point, tri_ind, tri_tree, nothing, nothing, nothing, material)
#     end
#     function RawMeshCache(point::Vector{SVector{3,Float64}}, tri_ind::Vector{SVector{3,Int64}}, tet_ind::Vector{SVector{4,Int64}}, strain::Vector{Float64}, material::MaterialProperties)
#         tri_tree = triTetMeshToTreeAABB(point, tri_ind)
#         tet_tree = triTetMeshToTreeAABB(point, tet_ind)
#         return new(point, tri_ind, tri_tree, tet_ind, tet_tree, strain, material)
#     end
# end
#
# function asHomogenousMesh(meshCache::RawMeshCache)
#     vec_Face = Face{3, Int32}.(meshCache.tri_ind)  # TODO: consider changing this to Int32
#     vec_Point = Point{3, Float32}.(meshCache.point)  # TODO: consider changing this to Float32
#     return HomogenousMesh(vec_Point, vec_Face)
# end
#
# @RigidBodyDynamics.indextype MeshID
# const MeshDict{V} = RigidBodyDynamics.IndexDict{MeshID, Base.OneTo{MeshID}, V}
# const MeshCacheDict{V} = RigidBodyDynamics.CacheIndexDict{MeshID, Base.OneTo{MeshID}, V}
# Base.@propagate_inbounds Base.getindex(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, key::RigidBody) = d[BodyID(key)]
# Base.@propagate_inbounds Base.setindex!(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, value, key::RigidBody) = d[BodyID(key)] = value
#
# struct MeshCache
#     name::String
#     BodyID::Union{Nothing,BodyID}
#     FrameID::CartesianFrame3D
#     raw::RawMeshCache
#     MeshCache(name::String, raw::RawMeshCache, body::RigidBody{Float64}) = new(name, BodyID(body), default_frame(body), raw)
# end
#
# new_mesh_cache() = MeshCacheDict{MeshCache}(Base.OneTo(MeshID(0)))
# function remake_mesh_cache(mesh_cache_old::MeshCacheDict{MeshCache}, name::String, raw::RawMeshCache, body::RigidBody{Float64})
#     mesh_n_old = length(mesh_cache_old)
#     mesh_n_new = mesh_n_old + 1
#     mesh_ids_old = Base.OneTo(MeshID(mesh_n_old))
#     mesh_ids_new = Base.OneTo(MeshID(mesh_n_new))
#     mesh_cache_new = MeshCacheDict{MeshCache}(mesh_ids_new)
#     for id = mesh_ids_old
#         mesh_cache_new[id] = mesh_cache_old[id]
#     end
#     mesh_cache_new[MeshID(mesh_n_new)] = MeshCache(name, raw, body)
#     return mesh_cache_new
# end
