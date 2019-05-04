
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

struct InertiaProperties{T<:Union{Tri,Tet}}
    d::Union{Nothing,Float64}  # if volume mesh is known thickness isn't needed to calculate inertia
    rho::Float64  # rho is always needed to calculate inertia
    function InertiaProperties(rho::Float64; d::Union{Nothing,Float64}=nothing)
        isa(d, Nothing) || (0.001 <= d <= 0.1) || error("thickness in unexpected range.")
        (50.0 <= rho) || error("rho in unexpected range.")
        i_prop_type = ifelse(d==nothing, Tet, Tri)
        return new{i_prop_type}(d, rho)
    end
end

struct MeshCache{T1,T2}
    name::String
    BodyID::Union{Nothing,BodyID}
    FrameID::CartesianFrame3D
    mesh::eMesh{T1,T2}
    tree::bin_BB_Tree{OBB}
    c_prop::Union{Nothing,ContactProperties}
    function MeshCache(name::String, e_mesh::eMesh{T1,T2}, tree::bin_BB_Tree, body::RigidBody{Float64},
        c_prop::Union{Nothing,ContactProperties}) where {T1,T2}

        return new{T1,T2}(name, BodyID(body), default_frame(body), e_mesh, tree, c_prop)
    end
end

get_tree(m::MeshCache) = m.tree

@inline get_c_prop(m::MeshCache) = m.c_prop
@inline Binary_BB_Trees.get_point(m::MeshCache) = m.mesh.point
@inline get_Ē(m::MeshCache) = get_c_prop(m).Ē
@inline get_ind_tri(m::MeshCache) = m.mesh.tri
@inline get_ind_tet(m::MeshCache) = m.mesh.tet
@inline get_ϵ(m::MeshCache) = m.mesh.ϵ

@RigidBodyDynamics.indextype MeshID
const MeshDict{V} = RigidBodyDynamics.IndexDict{MeshID, Base.OneTo{MeshID}, V}
const MeshCacheDict{V} = RigidBodyDynamics.CacheIndexDict{MeshID, Base.OneTo{MeshID}, V}
Base.@propagate_inbounds Base.getindex(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, key::RigidBody) = d[BodyID(key)]
Base.@propagate_inbounds Base.setindex!(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, value, key::RigidBody) = d[BodyID(key)] = value
