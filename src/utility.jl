
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

# findMesh(mech_scen::MechanismScenario, name::String) = findMesh(mech_scen.MeshCache, name)
# findMesh(ts::TempContactStruct, name::String) = findMesh(ts.MeshCache, name)
# findMesh(ts::MeshCacheDict{MeshCache}, name::String) = ts[findmesh(ts, name)]

function Radau_for_MechanismScenario(m::MechanismScenario{NX,NQ,Dual{Type_Tag,Float64,NC}}) where {NX,NQ,Type_Tag,NC}
    return makeRadauIntegrator(m, NX, 1.0e-16, 2, NC)
end

as_static_vector(f::Wrench{T}) where {T} = vcat(angular(f), linear(f))
as_static_vector(f::Twist{T}) where {T} = vcat(angular(f), linear(f))

@inline get_bristle_d0(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID) where {N,T} = segments(tm.s)[bristle_id]
@inline get_bristle_d1(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID) where {N,T} = segments(tm.ṡ)[bristle_id]
