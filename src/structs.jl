
struct MeshInertiaInfo
    tensor_I::SMatrix{3,3,Float64,9}
    com::SVector{3,Float64}
    mass::Float64
    mesh_vol::Float64
end

struct ContactProperties
    Ē::Float64
    μ::Float64
    d⁻¹::Float64
    χ::Float64
    function ContactProperties(;Ē::Float64, μ::Float64, d::Float64, χ::Float64)
        (1.0e4 <= Ē <= 1.0e8) || error("E_effective in unexpected range.")
        (0.0 <= μ <= 3.0) || error("mu in unexpected range.")
        (0.001 <= d <= 1.0) || error("thickness in unexpected range.")
        d⁻¹ = 1 / d
        (0.001 <= χ <= 5.0) || error("hc_velocity_damping in unexpected range.")
        return new(Ē, μ, d⁻¹, χ)
    end
end

calculateExtrensicCompliance(mat::ContactProperties) = 1 / (mat.Ē * mat.d⁻¹)

struct InertiaProperties
    d::Union{Nothing,Float64}  # if volume mesh is known thickness isn't needed to calculate inertia
    rho::Float64  # rho is always needed to calculate inertia
    function InertiaProperties(;rho::Union{Nothing,Float64}=nothing, d::Union{Nothing,Float64}=nothing)
        (rho == nothing) && error("rho is required")
        if d isa Float64
            (0.001 <= d < 0.1) || error("thickness in unexpected range.")
        end
        (50.0 <= rho < 2000.0) || error("rho in unexpected range.")
        return new(d, rho)
    end
end
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
    function TetMesh(point::Vector{SVector{3,Float64}}, tet_ind::Vector{SVector{4,Int64}}, ϵ::Vector{Float64},
            c_prop::ContactProperties)

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

    function MeshCache(name::String, h_mesh::HomogenousMesh, body::RigidBody{Float64},
            inertia_prop::Union{Nothing, InertiaProperties}=nothing)

        tri_simp_tree = SimplexTree(h_mesh)
        point = get_h_mesh_vertices(h_mesh)
        return new(point, name, BodyID(body), default_frame(body), inertia_prop, tri_simp_tree, nothing)
    end
end

@inline get_tri_mesh(m::MeshCache) = m.tri.ind
@inline get_tet_mesh(m::MeshCache) = m.tet.tet.ind
is_compliant(m::MeshCache) = (m.tet != nothing)

calc_mutual_μ(a::Nothing, b::Nothing) = error("both materials are nothing")
calc_mutual_μ(a::TetMesh, b::Nothing) = a.c_prop.μ
calc_mutual_μ(a::Nothing, b::TetMesh) = calc_mutual_μ(b, a)
calc_mutual_μ(a::TetMesh, b::TetMesh) = sqrt(a.c_prop.μ * b.c_prop.μ)
calc_mutual_μ(a::MeshCache, b::MeshCache) = calc_mutual_μ(a.tet, b.tet)

get_Ē(m::MeshCache) = m.tet.c_prop.Ē
