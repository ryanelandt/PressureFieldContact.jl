
function fill_with_nothing!(a)  # TODO: find more elegant way to do this
    for k = keys(a)
        a[k] = nothing
    end
    return nothing
end

@inline mat_mul_SA_bug_circumvent(A::SMatrix{4,4,Float64,16}, B::SMatrix{4,4,Float64,16}) = A * B
@inline mat_mul_SA_bug_circumvent(A::SMatrix{4,4,T,16}, B::SMatrix{4,4,T,16}) where {T} = StaticArrays.mul_loop(Size(A), Size(B), A, B)
@inline mat_mul_SA_bug_circumvent(A::SMatrix{4,4,Float64,16}, B::SMatrix{4,4,T,16}) where {T} = StaticArrays.mul_loop(Size(A), Size(B), A, B)
@inline mat_mul_SA_bug_circumvent(A::SMatrix{4,4,T,16}, B::SMatrix{4,4,Float64,16}) where {T} = StaticArrays.mul_loop(Size(A), Size(B), A, B)

function add_h_mesh_color(h_mesh::HomogenousMesh; color::Union{Nothing, RGBA{Float32}}=nothing)
    return HomogenousMesh(vertices=h_mesh.vertices, faces=h_mesh.faces, color=color)
end

function Radau_for_MechanismScenario(m::MechanismScenario{NX,NQ,Dual{Type_Tag,Float64,NC}}) where {NX,NQ,Type_Tag,NC}
    return makeRadauIntegrator(m, NX, 1.0e-16, 2, NC)
end

as_static_vector(f::Wrench{T}) where {T} = vcat(angular(f), linear(f))
as_static_vector(f::Twist{T}) where {T} = vcat(angular(f), linear(f))

function get_bristle_d0(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID) where {N,T}
    s = segments(tm.s)[bristle_id]
    return SVector{6,T}(s[1], s[2], s[3], s[4], s[5], s[6])
end
@inline get_bristle_d1(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID) where {N,T} = segments(tm.ṡ)[bristle_id]

function find_mesh_id(ts::MeshCacheDict{MeshCache}, name::String)  # TODO: make this function more elegant
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

function find_mesh_id(ts::MeshCacheDict{MeshCache}, mc::MeshCache)  # TODO: make this function more elegant
    id = MeshID(-9999)
    for k = keys(ts)
        if ts[k] == mc
            (id == MeshID(-9999)) || error("multiple")
            id = k
        end
    end
    (id == MeshID(-9999)) && error("no mesh found by name: $(mc.name)")
    return id
end

find_mesh_id(ts::TempContactStruct, input_2) = find_mesh_id(ts.MeshCache, input_2)

find_mesh(ts::MeshCacheDict{MeshCache}, name::String) = ts[find_mesh_id(ts, name)]
find_mesh(ts::Union{MechanismScenario,TempContactStruct}, name::String) = find_mesh(ts.MeshCache, name)
