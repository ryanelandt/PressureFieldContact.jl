struct SimplexTree{N}
    tree::bin_BB_Tree{AABB}
    ind::Vector{SVector{N,Int64}}
    function SimplexTree(h_mesh::HomogenousMesh)
        point = get_h_mesh_vertices(h_mesh)
        ind = get_h_mesh_faces(h_mesh)
        if length(ind) == 1
            tree = bin_BB_Tree{AABB}(1, svSvToAABB(point[ind[1]]))
        else
            tree = triTetMeshToTreeAABB(point, ind)
        end
        return new{3}(tree, ind)
    end
    function SimplexTree(point::Vector{SVector{3,Float64}}, ind::Vector{SVector{4,Int64}})
        tree = triTetMeshToTreeAABB(point, ind)
        return new{4}(tree, ind)
    end
end

struct TetMesh
    tet::SimplexTree{4}
    ϵ::Vector{Float64}
    c_prop::ContactProperties
    function TetMesh(point::Vector{SVector{3,Float64}}, tet_ind::Vector{SVector{4,Int64}}, ϵ::Vector{Float64}, c_prop::ContactProperties)
        tet_simp_tree = SimplexTree(point, tet_ind)
        return new(tet_simp_tree, ϵ, c_prop)
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

    function MeshCache(name::String, h_mesh::HomogenousMesh, tet_mesh::TetMesh, body::RigidBody{Float64},
            inertia_prop::Union{Nothing, InertiaProperties}=nothing)

        point, tri_ind = extract_HomogenousMesh_face_vertices(h_mesh)
        tri_simp_tree = SimplexTree(h_mesh)
        return new(point, name, BodyID(body), default_frame(body), inertia_prop, tri_simp_tree, tet_mesh)
    end

    function MeshCache(name::String, h_mesh::HomogenousMesh, body::RigidBody{Float64}, inertia_prop::Union{Nothing, InertiaProperties}=nothing)
        tri_simp_tree = SimplexTree(h_mesh)
        point = get_h_mesh_vertices(h_mesh)
        return new(point, name, BodyID(body), default_frame(body), inertia_prop, tri_simp_tree, nothing)
    end
end

function asHomogenousMesh(meshCache::MeshCache; color::Union{Nothing, RGBA{Float32}}=nothing)
    vec_Face = Face{3, Int32}.(meshCache.tri.ind)
    vec_Point = Point{3, Float32}.(meshCache.point)
    return HomogenousMesh(vertices=vec_Point, faces=vec_Face, color=color)
end

@RigidBodyDynamics.indextype MeshID
const MeshDict{V} = RigidBodyDynamics.IndexDict{MeshID, Base.OneTo{MeshID}, V}
const MeshCacheDict{V} = RigidBodyDynamics.CacheIndexDict{MeshID, Base.OneTo{MeshID}, V}
Base.@propagate_inbounds Base.getindex(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, key::RigidBody) = d[BodyID(key)]
Base.@propagate_inbounds Base.setindex!(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, value, key::RigidBody) = d[BodyID(key)] = value
