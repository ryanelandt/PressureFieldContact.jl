struct RawMeshCache
    point::Vector{SVector{3,Float64}}
    tri_ind::Vector{SVector{3,Int64}}
    tri_tree::bin_BB_Tree{AABB}
    tet_ind::Union{Nothing,Vector{SVector{4,Int64}}}
    tet_tree::Union{Nothing,bin_BB_Tree{AABB}}
    strain::Union{Nothing,Vector{Float64}}
    material::MaterialProperties
    function RawMeshCache(point::Vector{SVector{3,Float64}}, tri_ind::Vector{SVector{3,Int64}}, material::MaterialProperties=MaterialProperties())
        tri_tree = triTetMeshToTreeAABB(point, tri_ind)
        return new(point, tri_ind, tri_tree, nothing, nothing, nothing, material)
    end
    function RawMeshCache(point::Vector{SVector{3,Float64}}, tri_ind::Vector{SVector{3,Int64}}, tet_ind::Vector{SVector{4,Int64}}, strain::Vector{Float64}, material::MaterialProperties)
        tri_tree = triTetMeshToTreeAABB(point, tri_ind)
        tet_tree = triTetMeshToTreeAABB(point, tet_ind)
        return new(point, tri_ind, tri_tree, tet_ind, tet_tree, strain, material)
    end
end

@RigidBodyDynamics.indextype MeshID
const MeshDict{V} = RigidBodyDynamics.IndexDict{MeshID, Base.OneTo{MeshID}, V}
const MeshCacheDict{V} = RigidBodyDynamics.CacheIndexDict{MeshID, Base.OneTo{MeshID}, V}
Base.@propagate_inbounds Base.getindex(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, key::RigidBody) = d[BodyID(key)]
Base.@propagate_inbounds Base.setindex!(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, value, key::RigidBody) = d[BodyID(key)] = value

struct MeshCache
    name::String
    BodyID::Union{Nothing,BodyID}
    raw::RawMeshCache
    MeshCache(name::String, raw::RawMeshCache, body::RigidBody{Float64}) = new(name, BodyID(body), raw)
end
