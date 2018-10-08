### RigidBodyDynamics ###
function RigidBodyDynamics.principal_value!(mech_scen::MechanismScenario, x::Vector{Float64})
    copyto!(mech_scen.float, x)
    principal_value!(mech_scen.float.state)
    copyto!(x, mech_scen.float)
    return nothing
end

### NumericalTricks ###
NumericalTricks.safe_normalize(p::FreeVector3D{SVector{3,T}}) where {T} = FreeVector3D(p.frame, safe_normalize(p.v))
function NumericalTricks.vec_sub_vec_proj(vec_in::FreeVector3D{SVector{3, TD}}, n̂::FreeVector3D{SVector{3,TD}}) where {TD}
    @framecheck(vec_in.frame, n̂.frame)
    return FreeVector3D(vec_in.frame, vec_sub_vec_proj(vec_in.v, n̂.v))
end
function NumericalTricks.soft_clamp(fv::FreeVector3D{SVector{3,T}}, bound::FreeVector3D{SVector{3,T}}) where {T}  # TODO soften this
    @framecheck(fv.frame, bound.frame)
    abs_bound = abs.(bound.v)
    fv_clamp = soft_clamp.(fv.v, abs_bound)
    return FreeVector3D(bound.frame, fv_clamp)
end
@inline NumericalTricks.norm_squared(fv::FreeVector3D{SVector{3,T}}) where {T} = norm_squared(fv.v)

### Tri_Tet_Intersections ###
Tri_Tet_Intersections.area(p1::Point3D{T}, p2::Point3D{T}, p3::Point3D{T}) where {T} = area(p1.v, p2.v, p3.v)
Tri_Tet_Intersections.getTop(m::MatrixTransform{4,N2,T,N3}) where {N2,T,N3} = MatrixTransform(m.from, m.to, getTop(m.mat))
function Tri_Tet_Intersections.add!(c::ClippedPolygon{4,T}, p::Point4D{SVector{4,T}}) where {T}
    @framecheck(c.frame, p.frame)
    add!(c, p.v)
    return nothing
end
function Tri_Tet_Intersections.tet_clip_poly_to_cartesian!(poly_3D::ClippedPolygon{3,T}, poly_4D_1::ClippedPolygon{4,T}, A_w_zeta_top::MatrixTransform) where{T}
    @framecheck(poly_4D_1.frame, A_w_zeta_top.from)
    @framecheck(poly_3D.frame, A_w_zeta_top.to)
    tet_clip_poly_to_cartesian!(poly_3D, poly_4D_1, A_w_zeta_top.mat)
    return nothing
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
