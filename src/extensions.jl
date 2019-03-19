### RigidBodyDynamics ###
function RigidBodyDynamics.principal_value!(mech_scen::MechanismScenario, x::Vector{Float64})
    copyto!(mech_scen.float, x)
    principal_value!(mech_scen.float.state)
    copyto!(x, mech_scen.float)
    return nothing
end

function RigidBodyDynamics.set_configuration!(mech_scen::MechanismScenario, joint::Joint, config)
    return set_configuration!(mech_scen.float.state, joint, config)
end

### NumericalTricks ###
NumericalTricks.safe_normalize(p::FreeVector3D{SVector{3,T}}) where {T} = FreeVector3D(p.frame, safe_normalize(p.v))
function NumericalTricks.vec_sub_vec_proj(vec_in::FreeVector3D{SVector{3, TD}}, n̂::FreeVector3D{SVector{3,TD}}) where {TD}
    @framecheck(vec_in.frame, n̂.frame)
    return FreeVector3D(vec_in.frame, vec_sub_vec_proj(vec_in.v, n̂.v))
end
@inline NumericalTricks.norm²(fv::FreeVector3D{SVector{3,T}}) where {T} = norm²(fv.v)

### Tri_Tet_Intersections ###
Tri_Tet_Intersections.area(p1::Point3D{T}, p2::Point3D{T}, p3::Point3D{T}) where {T} = area(p1.v, p2.v, p3.v)
Tri_Tet_Intersections.getTop(m::MatrixTransform{4,N2,T,N3}) where {N2,T,N3} = MatrixTransform(m.from, m.to, getTop(m.mat))
Tri_Tet_Intersections.unPad(p::Point4D{SVector{4,T}}) where {T} = Point3D(p.frame, SVector{3,T}(p.v[1], p.v[2], p.v[3]))

### Binary_BB_Trees
Binary_BB_Trees.fit_tri_obb(eM::eMesh{Tri,T2}, k::Int64) where {T2} = fit_tri_obb(eM.point[eM.tri[k]])
function Binary_BB_Trees.fit_tet_obb(eM::eMesh{Tri,T2}, k::Int64) where {T2}
    i = eM.tet[k]
    return fit_tet_obb(eM.point[i], eM.ϵ[i])
end

### Base ###
function Base.copyto!(dest::TypedMechanismScenario{N,T}, src::AbstractVector{T}) where {N,T}
    nq = num_positions(dest.state)
    nv = num_velocities(dest.state)
    ns = length(dest.s)
    @boundscheck length(src) == nq + nv + ns || throw(DimensionMismatch())
    @inbounds copyto!(parent(dest.state.q), 1, src, 1, nq)
    @inbounds copyto!(parent(dest.state.v), 1, src, nq + 1, nv)
    @inbounds copyto!(dest.s, 1, src, nq + nv + 1, ns)
    setdirty!(dest.state)
    dest
end
function Base.copyto!(dest::AbstractVector{T}, src::TypedMechanismScenario{N,T}) where {N,T}
    nq = num_positions(src.state)
    nv = num_velocities(src.state)
    ns = length(src.s)
    length(dest) == nq + nv + ns || throw(DimensionMismatch())
    @inbounds copyto!(dest, 1, src.state.q, 1, nq)
    @inbounds copyto!(dest, nq + 1, src.state.v, 1, nv)
    @inbounds copyto!(dest, nq + nv + 1, src.s, 1, ns)
    dest
end
function Base.copyto!(ẋ::AbstractVector, src::TypedMechanismScenario{N,T}, result::DynamicsResult) where {N,T}
    nq = length(result.q̇)
    nv = length(result.v̇)
    ns = length(src.ṡ)
    @boundscheck length(ẋ) == nq + nv + ns || throw(DimensionMismatch())
    @inbounds copyto!(ẋ, 1, result.q̇, 1, nq)
    @inbounds copyto!(ẋ, nq + 1, result.v̇, 1, nv)
    @inbounds copyto!(ẋ, nq + nv + 1, src.ṡ, 1, ns)
    ẋ
end
