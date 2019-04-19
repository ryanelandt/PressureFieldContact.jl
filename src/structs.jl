
struct MeshInertiaInfo
    tensor_I::SMatrix{3,3,Float64,9}
    com::SVector{3,Float64}
    mass::Float64
    mesh_vol::Float64
end

struct ContactProperties
    Ē::Float64
    function ContactProperties(;Ē::Float64)
        (1.0e4 <= Ē <= 3.0e11) || error("E_effective in unexpected range.")
        return new(Ē)
    end
end

struct eTree{T1<:Union{Nothing,Tri},T2<:Union{Nothing,Tet}}
    tri::Union{Nothing,bin_BB_Tree}
    tet::Union{Nothing,bin_BB_Tree}
    c_prop::Union{Nothing,ContactProperties}
    function eTree(e_mesh::eMesh{Tri,Tet}, c_prop::ContactProperties)
        tri = triTetMeshToTreeAABB(e_mesh.point, e_mesh.tri)
        tet = triTetMeshToTreeAABB(e_mesh.point, e_mesh.tet)
        return new{Tri,Tet}(tri, tet, c_prop)
    end
    function eTree(e_mesh::eMesh{Tri,Nothing}, c_prop::Nothing=nothing)
        tri = triTetMeshToTreeAABB(e_mesh.point, e_mesh.tri)
        return new{Tri,Nothing}(tri, nothing, nothing)
    end
    function eTree(e_mesh::eMesh{Nothing,Tet}, c_prop::ContactProperties)
        tet = triTetMeshToTreeAABB(e_mesh.point, e_mesh.tet)
        return new{Nothing,Tet}(nothing, tet, c_prop)
    end
    ###
    eTree(tri::Nothing,     tet::bin_BB_Tree, c_prop::ContactProperties) = new{Nothing,Tet}(tri, tet, c_prop)
    eTree(tri::bin_BB_Tree, tet::Nothing,     c_prop::Nothing)           = new{Tri,Nothing}(tri, tet, c_prop)
    eTree(tri::bin_BB_Tree, tet::bin_BB_Tree, c_prop::ContactProperties) = new{Tri,Tet}(tri, tet, c_prop)
end

struct InertiaProperties
    d::Union{Nothing,Float64}  # if volume mesh is known thickness isn't needed to calculate inertia
    rho::Float64  # rho is always needed to calculate inertia
    function InertiaProperties(rho::Float64; d::Union{Nothing,Float64}=nothing)
        isa(d, Nothing) || (0.001 <= d <= 0.1) || error("thickness in unexpected range.")
        (50.0 <= rho) || error("rho in unexpected range.")
        return new(d, rho)
    end
end

struct MeshCache{T1,T2}
    name::String
    BodyID::Union{Nothing,BodyID}
    FrameID::CartesianFrame3D
    mesh::eMesh{T1,T2}
    tree::eTree{T1,T2}
    function MeshCache(name::String, e_mesh::eMesh{T1,T2}, e_tree::eTree{T1,T2}, body::RigidBody{Float64}) where {T1,T2}
        return new{T1,T2}(name, BodyID(body), default_frame(body), e_mesh, e_tree)
    end
end

get_tree(m::MeshCache{Tri,Nothing}) = m.tree.tri
get_tree(m::MeshCache{Nothing,Tet}) = m.tree.tet

# @inline get_tree_tet(m::MeshCache) = m.tree.tet
# @inline get_tree_tri(m::MeshCache) = m.tree.tri

@inline get_c_prop(m::MeshCache) = m.tree.c_prop
@inline Binary_BB_Trees.get_point(m::MeshCache) = m.mesh.point
@inline get_ind_tri(m::MeshCache) = m.mesh.tri
@inline get_ind_tet(m::MeshCache) = m.mesh.tet
@inline get_ϵ(m::MeshCache) = m.mesh.ϵ

is_compliant(m::MeshCache{T1,Tet}) where {T1} = true
is_compliant(m::MeshCache{T1,Nothing}) where {T1} = false

is_rigid(m::MeshCache{Tri,T2}) where {T2} = true
is_rigid(m::MeshCache{Nothing,T2}) where {T2} = false

get_Ē(m::MeshCache) = get_c_prop(m).Ē

@RigidBodyDynamics.indextype MeshID
const MeshDict{V} = RigidBodyDynamics.IndexDict{MeshID, Base.OneTo{MeshID}, V}
const MeshCacheDict{V} = RigidBodyDynamics.CacheIndexDict{MeshID, Base.OneTo{MeshID}, V}
Base.@propagate_inbounds Base.getindex(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, key::RigidBody) = d[BodyID(key)]
Base.@propagate_inbounds Base.setindex!(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, value, key::RigidBody) = d[BodyID(key)] = value
