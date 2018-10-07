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

function asHomogenousMesh(meshCache::MeshCache)
    vec_Face = Face{3, Int32}.(meshCache.tri.ind)
    vec_Point = Point{3, Float32}.(meshCache.point)
    return HomogenousMesh(vec_Point, vec_Face)
end

asHomogenousMesh(meshCache::MeshCache, color) = asHomogenousMesh(meshCache, RGBA{Float32}(color...))
function asHomogenousMesh(meshCache::MeshCache, color::RGBA{Float32})
    vec_Face = Face{3, Int32}.(meshCache.tri.ind)
    vec_Point = Point{3, Float32}.(meshCache.point)
    return HomogenousMesh(vec_Point, vec_Face, color)
end

@RigidBodyDynamics.indextype MeshID
const MeshDict{V} = RigidBodyDynamics.IndexDict{MeshID, Base.OneTo{MeshID}, V}
const MeshCacheDict{V} = RigidBodyDynamics.CacheIndexDict{MeshID, Base.OneTo{MeshID}, V}
Base.@propagate_inbounds Base.getindex(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, key::RigidBody) = d[BodyID(key)]
Base.@propagate_inbounds Base.setindex!(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, value, key::RigidBody) = d[BodyID(key)] = value
