
function calcXd!(xx::AbstractVector{T}, x::AbstractVector{T}, m::MechanismScenario{NX,NQ,T}, t::Float64=0.0) where {NX,NQ,T}
    return calcXd!(xx, x, m, m.dual, t)
end

function calcXd!(xx::AbstractVector{Float64}, x::AbstractVector{Float64}, m::MechanismScenario{NX,NQ,T}, t::Float64=0.0) where {NX,NQ,T}
    return calcXd!(xx, x, m, m.float, t)
end

```
Conventions:
n̂ refers to the contact surface normal that points into body B
v_cart refers to + velocity of B - the velocity of A
the wrench is the wrench applied TO body A
```
function calcXd!(xx::AbstractVector{T1}, x::AbstractVector{T1}, m::MechanismScenario{NX,NQ,T2},
        tm::TypedMechanismScenario{NQ,T1}, t::Float64=0.0) where {NX,NQ,T1,T2}

    state = tm.state
    copyto!(tm, x)
    H = tm.result.massmatrix
    mass_matrix!(H, state)
    dynamics_bias!(tm.result, state)
    configuration_derivative!(tm.result.q̇, state)
    forceAllElasticIntersections!(m, tm)

    (m.continuous_controller == nothing) || m.continuous_controller(tm, t)

	sum_all_forces!(m, tm)

	chol_fact = LinearAlgebra.cholesky!(H)
	ldiv!(tm.result.v̇.parent, chol_fact, tm.rhs)

    copyto!(xx, tm, tm.result)
    return nothing
end

function sum_all_forces!(m::MechanismScenario{NX,NQ,T2}, tm::TypedMechanismScenario{NQ,Float64}) where {NX,NQ,T2}
	BLAS.blascopy!(tm.nv, tm.f_generalized, 1, tm.rhs, 1)
	BLAS.axpy!(-1.0, tm.result.dynamicsbias.parent, tm.rhs)
	BLAS.axpy!(+1.0, tm.τ_ext, tm.rhs)
	BLAS.axpy!(+1.0, m.τ_ext, tm.rhs)
end

function sum_all_forces!(m::MechanismScenario{NX,NQ,T2}, tm::TypedMechanismScenario{NQ,T1}) where {NX,NQ,T1,T2}
	tm.rhs  .= tm.f_generalized
	tm.rhs .-= tm.result.dynamicsbias.parent
	tm.rhs .+= tm.τ_ext
	tm.rhs .+= m.τ_ext
end

function calcXd(x::AbstractVector{T1}, m::MechanismScenario{NX,NQ,T2}, t::Float64=0.0) where {T1,NX,NQ,T2}
    xx = deepcopy(x)
    calcXd!(xx, x, m)
    return xx
end

function forceAllElasticIntersections!(m::MechanismScenario{NX,NQ,T1}, tm::TypedMechanismScenario{NQ,T2}) where {NX,NQ,T1,T2}
    refreshJacobians!(m, tm)
    tm.f_generalized .= zero(T2)
    for k = 1:length(m.ContactInstructions)
        c_ins = m.ContactInstructions[k]
        force_single_elastic_intersection!(m, tm, c_ins)
    end
    return nothing
end

function force_single_elastic_intersection!(m::MechanismScenario{NX,NQ,T1}, tm::TypedMechanismScenario{NQ,T2},
    c_ins::ContactInstructions) where {NX,NQ,T1,T2}

    calcTriTetIntersections!(m, c_ins)
    if (0 != length(m.TT_Cache))  # yes intersections
        refreshBodyBodyCache!(m, tm, c_ins)
        integrate_over_logic!(tm.bodyBodyCache, m.TT_Cache)
        if !isempty(tm.bodyBodyCache.TractionCache)
            wrench = yes_contact!(c_ins.FrictionModel, tm, c_ins)
            addGeneralizedForcesThirdLaw!(wrench, tm, c_ins)
            return wrench
        end
    end
    no_contact!(c_ins.FrictionModel, tm, c_ins)
end

function refreshJacobians!(m::MechanismScenario{NX,NQ,T1}, tm::TypedMechanismScenario{NQ,T2}) where {NX,NQ,T1,T2}
    for k = m.body_ids
        path_k = m.path[k]
        (path_k != nothing) && geometric_jacobian!(tm.GeometricJacobian[k], tm.state, path_k)
    end
    return nothing
end

function calcTriTetIntersections!(m::MechanismScenario, c_ins::ContactInstructions) # where {N,T}
    b = m.float.bodyBodyCache  # this can be float because intersection is assumed to not depend on partials
    refreshBodyBodyTransform!(m, m.float, c_ins)  # TODO: isn't b.x_r¹_rʷ already calculated?
    x_r¹_r² = inv(b.x_rʷ_r¹) * b.x_rʷ_r²
    update_TT_Cache!(m.TT_Cache, translation(x_r¹_r²), rotation(x_r¹_r²))
    tree_2 = get_tree_tet(b.mesh_2)
    tree_1 = ifelse(c_ins.mutual_compliance, get_tree_tet, get_tree_tri)(b.mesh_1)
    tree_tree_intersect(m.TT_Cache, tree_1, tree_2)
    return nothing
end

function refreshBodyBodyTransform!(m::MechanismScenario, tm::TypedMechanismScenario{N,T},
        c_ins::ContactInstructions) where {N,T}

    b = tm.bodyBodyCache
    b.mesh_1 = m.MeshCache[c_ins.id_1]
    b.mesh_2 = m.MeshCache[c_ins.id_2]
    b.x_rʷ_r¹ = transform_to_root(tm.state, b.mesh_1.BodyID)  # TODO: add safe=false
    b.x_rʷ_r² = transform_to_root(tm.state, b.mesh_2.BodyID)  # TODO: add safe=false
    b.x_r²_rʷ = inv(b.x_rʷ_r²)
    b.x_r¹_rʷ = inv(b.x_rʷ_r¹)
    return nothing
end

function refreshBodyBodyCache!(m::MechanismScenario, tm::TypedMechanismScenario{N,T},
        c_ins::ContactInstructions) where {N,T}

    b = tm.bodyBodyCache
    empty!(b.TractionCache)
    refreshBodyBodyTransform!(m, tm, c_ins)

    # Rates
    twist_w_r¹ = twist_wrt_world(tm.state, b.mesh_1.BodyID)
    twist_w_r² = twist_wrt_world(tm.state, b.mesh_2.BodyID)
    b.twist_r¹_r² = -twist_w_r² + twist_w_r¹  # velocity of tri wrt tet expressed in world

    b.μ = c_ins.μ
    b.χ = c_ins.χ

    c_prop = get_c_prop(b.mesh_2)
    b.Ē = c_prop.Ē
    b.d⁻¹ = c_prop.d⁻¹
    return nothing
end

function integrate_over_logic!(b::TypedElasticBodyBodyCache{N,T}, ttCache::TT_Cache) where {N,T}
    mesh_1 = b.mesh_1
    mesh_2 = b.mesh_2
    x_rʷ_r¹ = b.x_rʷ_r¹
    x_rʷ_r² = b.x_rʷ_r²
    x_r¹_rʷ = b.x_r¹_rʷ
    x_r²_rʷ = b.x_r²_rʷ
    if is_compliant(mesh_1)  # this is not right
        integrate_over_volume_volume_all!(mesh_1, mesh_2, x_rʷ_r¹, x_rʷ_r², x_r¹_rʷ, x_r²_rʷ, b, ttCache)
    else
        integrate_over_surface_volume_all!(mesh_1, mesh_2, x_rʷ_r¹, x_rʷ_r², x_r¹_rʷ, x_r²_rʷ, b, ttCache)
    end
end

function integrate_over_volume_volume_all!(mesh_1::MeshCache, mesh_2::MeshCache, x_rʷ_r¹::Transform3D{T},
        x_rʷ_r²::Transform3D{T}, x_r¹_rʷ::Transform3D{T}, x_r²_rʷ::Transform3D{T}, b::TypedElasticBodyBodyCache{N,T},
        ttCache::TT_Cache) where {N,T}

    for k = 1:length(ttCache.vc)
        i_1, i_2 = ttCache.vc[k]
        integrate_over_volume_volume!(i_1, i_2, mesh_1, mesh_2, x_rʷ_r¹, x_rʷ_r², x_r¹_rʷ, x_r²_rʷ, b)
    end
end

function triangle_vertices(i_tri::Int64, m::MeshCache)
    ind_vert = get_ind_tri(m)[i_tri]
    return get_point(m)[ind_vert]
end

function tetrahedron_vertices_ϵ(i_tet::Int64, m::MeshCache)
    ind_vert = get_ind_tet(m)[i_tet]
    ϵ = get_ϵ(m)[ind_vert]
    ϵ = SMatrix{1,4,Float64,4}(ϵ)
    cart_vert = get_point(m)[ind_vert]
    return cart_vert, ϵ
end

function calc_ζ_transforms(frame_ζ::CartesianFrame3D, frame_b ::CartesianFrame3D, p_tet, x_r¹_rʷ, x_rʷ_r¹)
    x_r¹_ζ = MatrixTransform(frame_ζ, frame_b, asMatOnePad(p_tet))
    x_ζ_r¹ = inv(x_r¹_ζ)
    x_rʷ_ζ = x_rʷ_r¹ * x_r¹_ζ
    x_ζ_rʷ = x_ζ_r¹ * x_r¹_rʷ  # NOTE: inv(x_r_ζ) is **always** Float64
    return x_rʷ_ζ, x_ζ_rʷ, x_ζ_r¹
end

find_plane_tet(E::Float64, ϵ::SMatrix{1,4,Float64,4}, X_r_w) = (E * ϵ) * X_r_w

function integrate_over_volume_volume!(i_1::Int64, i_2::Int64, mesh_1::MeshCache, mesh_2::MeshCache,
        x_rʷ_r¹::Transform3D{T}, x_rʷ_r²::Transform3D{T}, x_r¹_rʷ::Transform3D{T}, x_r²_rʷ::Transform3D{T},
        b::TypedElasticBodyBodyCache{N,T}) where {N,T}

    vert_1, ϵ¹ = tetrahedron_vertices_ϵ(i_1, mesh_1)
    vert_2, ϵ² = tetrahedron_vertices_ϵ(i_2, mesh_2)
    x_rʷ_ζ¹, x_ζ¹_rʷ, x_ζ¹_r¹ = calc_ζ_transforms(FRAME_ζ¹, mesh_1.FrameID, vert_1, x_r¹_rʷ, x_rʷ_r¹)
    x_rʷ_ζ², x_ζ²_rʷ, x_ζ²_r² = calc_ζ_transforms(FRAME_ζ², mesh_2.FrameID, vert_2, x_r²_rʷ, x_rʷ_r²)

    # TODO: make this better
    ϵ_plane_1_rʷ = find_plane_tet(get_Ē(mesh_1), ϵ¹, x_ζ¹_rʷ.mat)
    ϵ_plane_2_rʷ = find_plane_tet(get_Ē(mesh_2), ϵ², x_ζ²_rʷ.mat)
	ϵ_plane_rʷ = ϵ_plane_2_rʷ - ϵ_plane_1_rʷ  # normalize penetration extent is positive so this describes the plane
		# of the contact surface pointing towards mesh_2

    poly_rʷ = clip_plane_tet(ϵ_plane_rʷ, x_rʷ_ζ¹.mat)
    if 3 <= length(poly_rʷ)
        poly_ζ² = one_pad_then_mul(x_ζ²_rʷ.mat, poly_rʷ)
        poly_ζ² = zero_small_coordinates(poly_ζ²)  # This needs to be done to avoid a degenerate situation where the
            # plane lies exactly on the intersection of the faces of two tets. This situation happens EVERY time two
            # tet faces that lie on the surface intersect.
        poly_ζ² = clip_in_tet_coordinates(poly_ζ²)
        if 3 <= length(poly_ζ²)
			frame_world = b.frame_world
			n = unPad(ϵ_plane_rʷ)
			n̂ = unsafe_normalize(n)
            n̂ = FreeVector3D(frame_world, n̂)
            integrate_over_polygon_patch!(b, poly_ζ², frame_world, n̂, x_rʷ_ζ², x_ζ²_rʷ, ϵ², x_ζ²_r²)
        end
    end
end

function integrate_over_surface_volume_all!(mesh_1::MeshCache, mesh_2::MeshCache, x_rʷ_r¹::Transform3D{T},
        x_rʷ_r²::Transform3D{T}, x_r¹_rʷ::Transform3D{T}, x_r²_rʷ::Transform3D{T},
        b::TypedElasticBodyBodyCache{N,T}, ttCache::TT_Cache) where {N,T}

    for k = 1:length(ttCache.vc)
        i_1, i_2 = ttCache.vc[k]
        integrate_over_surface_volume!(i_1, i_2, mesh_1, mesh_2, x_rʷ_r¹, x_rʷ_r², x_r¹_rʷ, x_r²_rʷ, b)
    end
end

function integrate_over_surface_volume!(i_1::Int64, i_2::Int64, mesh_1::MeshCache, mesh_2::MeshCache,
        x_rʷ_r¹::Transform3D{T}, x_rʷ_r²::Transform3D{T}, x_r¹_rʷ::Transform3D{T}, x_r²_rʷ::Transform3D{T},
        b::TypedElasticBodyBodyCache{N,T}) where {N,T}

    vert_1 = triangle_vertices(i_1, mesh_1)
    vert_2, ϵ² = tetrahedron_vertices_ϵ(i_2, mesh_2)
    x_rʷ_ζ², x_ζ²_rʷ, x_ζ²_r² = calc_ζ_transforms(FRAME_ζ², mesh_2.FrameID, vert_2, x_r²_rʷ, x_rʷ_r²)
    x_ζ²_r¹ = x_ζ²_rʷ * x_rʷ_r¹

    n̂_r¹ = FreeVector3D(mesh_1.FrameID, triangleNormal(vert_1))  # REDACTED: pressure is applied opposite the triangle normal
    n̂_rʷ = x_rʷ_r¹ * n̂_r¹

    poly_rʷ = poly_eight(vert_1.data)
    poly_ζ² = one_pad_then_mul(x_ζ²_r¹.mat, poly_rʷ)
    poly_ζ² = clip_in_tet_coordinates(poly_ζ²)
    if 3 <= length(poly_ζ²)
        frame_world = b.frame_world
        integrate_over_polygon_patch!(b, poly_ζ², frame_world, n̂_rʷ, x_rʷ_ζ², x_ζ²_rʷ, ϵ², x_ζ²_r²)
    end
end

function integrate_over_polygon_patch!(b::TypedElasticBodyBodyCache{N,T}, poly_ζ²::poly_eight{4,T},
        frame_world::CartesianFrame3D, n̂::FreeVector3D{SVector{3,T}}, x_rʷ_ζ²::MatrixTransform{4,4,T,16},
        x_ζ²_rʷ::MatrixTransform{4,4,T,16}, ϵ::SMatrix{1,4,Float64,4}, x_ζ²_r²::MatrixTransform{4,4,Float64,16}) where {N,T}

    poly_rʷ = mul_then_un_pad(x_rʷ_ζ².mat, poly_ζ²)
    centroid_rʷ = Point3D(frame_world, centroid(poly_rʷ, n̂.v)[2])
    centroid_ζ² = x_ζ²_rʷ * centroid_rʷ
    N_vertices = length(poly_ζ²)
    ζ²_2 = getPoint(poly_ζ², FRAME_ζ², N_vertices)
    vert_rʷ_2 = getPoint(poly_rʷ, frame_world, N_vertices)
    for k = 1:N_vertices
        ζ²_1 = ζ²_2
        ζ²_2 = getPoint(poly_ζ², FRAME_ζ², k)
        x_ζ²_ϕ = hcat(ζ²_1, ζ²_2, centroid_ζ²)
        vert_rʷ_1 = vert_rʷ_2
        vert_rʷ_2 = getPoint(poly_rʷ, frame_world, k)
		# area_quad_k = area(vert_rʷ_1, vert_rʷ_2, centroid_rʷ, n̂)
		area_quad_k = triangle_area((vert_rʷ_1.v, vert_rʷ_2.v, centroid_rʷ.v), n̂.v)
		(-1.0e-6 <= area_quad_k) || error("area is negative!!! $(area_quad_k)")
        (0.0 < area_quad_k) && fillTractionCacheForTriangle!(b, area_quad_k, x_ζ²_ϕ, x_rʷ_ζ², n̂, ϵ, x_ζ²_r²)  # no point in adding intersection if area is zero
    end
    return nothing
end

# # TODO: create fillTractionCacheForTriangle! with a macro
function fillTractionCacheForTriangle!(b::TypedElasticBodyBodyCache{1,T}, area_quad_k::T,
        A_ζ_ϕ::MatrixTransform{4,3,T,12}, A_w_ζ::MatrixTransform{4,4,T,16}, n̂::FreeVector3D{SVector{3,T}},
        ϵ::SMatrix{1,4,Float64,4}, x_ζ²_r²::MatrixTransform{4,4,Float64,16}) where {T}

    r_cart_1, v_cart_t_1, dA_1, p_1 = fillTractionCacheInnerLoop!(1, b, A_ζ_ϕ, A_w_ζ, n̂, area_quad_k, ϵ, x_ζ²_r²)
    p = (p_1, )
    if 0.0 < sum(p)
        r_cart = (r_cart_1, )
        v_cart_t = (v_cart_t_1, )
        dA = (dA_1, )
        trac_cache = TractionCache(n̂, r_cart, v_cart_t, dA, p)
        addCacheItem!(b.TractionCache, trac_cache)
    end
    return nothing
end

function fillTractionCacheForTriangle!(b::TypedElasticBodyBodyCache{3,T}, area_quad_k::T,
        A_ζ_ϕ::MatrixTransform{4,3,T,12}, A_w_ζ::MatrixTransform{4,4,T,16}, n̂::FreeVector3D{SVector{3,T}},
        ϵ::SMatrix{1,4,Float64,4}, x_ζ²_r²::MatrixTransform{4,4,Float64,16}) where {T}

    r_cart_1, v_cart_t_1, dA_1, p_1 = fillTractionCacheInnerLoop!(1, b, A_ζ_ϕ, A_w_ζ, n̂, area_quad_k, ϵ, x_ζ²_r²)
    r_cart_2, v_cart_t_2, dA_2, p_2 = fillTractionCacheInnerLoop!(2, b, A_ζ_ϕ, A_w_ζ, n̂, area_quad_k, ϵ, x_ζ²_r²)
    r_cart_3, v_cart_t_3, dA_3, p_3 = fillTractionCacheInnerLoop!(3, b, A_ζ_ϕ, A_w_ζ, n̂, area_quad_k, ϵ, x_ζ²_r²)
    p = (p_1, p_2, p_3)
    if 0.0 < sum(p)
        dA = (dA_1, dA_2, dA_3)
        r_cart = (r_cart_1, r_cart_2, r_cart_3)
        v_cart_t = (v_cart_t_1, v_cart_t_2, v_cart_t_3)
        trac_cache = TractionCache(n̂, r_cart, v_cart_t, dA, p)
        addCacheItem!(b.TractionCache, trac_cache)
    end
    return nothing
end

function fillTractionCacheInnerLoop!(k::Int64, b::TypedElasticBodyBodyCache{N,T},
        A_ζ_ϕ::MatrixTransform{4,3,T,12}, A_w_ζ::MatrixTransform{4,4,T,16}, n̂::FreeVector3D{SVector{3,T}},
        area_quad_k::T, ϵ_tet::SMatrix{1,4,Float64,4}, x_ζ²_r²::MatrixTransform{4,4,Float64,16}) where {N,T}

    quad_point_ϕ = Point3D(FRAME_ϕ, b.quad.zeta[k])
    quad_point_ζ = A_ζ_ϕ * quad_point_ϕ
    p_cart_qp = unPad(A_w_ζ * quad_point_ζ)
    ϵ_quad = dot(ϵ_tet, quad_point_ζ.v)
    ṙ = point_velocity(b.twist_r¹_r², p_cart_qp)
    ṙ² = transform(ṙ, b.x_r²_rʷ)
    ϵϵ = dot(ϵ_tet, (x_ζ²_r² * ṙ²).v)
    damp_term = fastSoftPlus(1.0 + b.χ * ϵϵ)
    p = ϵ_quad * b.Ē
    p_hc = p * damp_term
    dA = b.quad.w[k] * area_quad_k
    return p_cart_qp, ṙ, dA, p_hc
end

function addGeneralizedForcesThirdLaw!(wrench::Wrench{T}, tm::TypedMechanismScenario{N,T}, cInfo::ContactInstructions) where {N,T}
    torque_third_law = tm.torque_third_law
    f_generalized = tm.f_generalized
    addGeneralizedForcesExternal!(f_generalized,  wrench, tm.GeometricJacobian[tm.bodyBodyCache.mesh_1.BodyID], torque_third_law)
    addGeneralizedForcesExternal!(f_generalized, -wrench, tm.GeometricJacobian[tm.bodyBodyCache.mesh_2.BodyID], torque_third_law)
end

function addGeneralizedForcesExternal!(f_generalized::Vector{T}, wrench::Wrench{T}, jac::Nothing,
        torque_third_law::Vector{T}) where {T}

    return nothing
end
function addGeneralizedForcesExternal!(f_generalized::Vector{T}, wrench::Wrench{T}, jac::GeometricJacobian{Matrix{T}},
        torque_third_law::Vector{T}) where {T}

    torque!(torque_third_law, jac, wrench)
    f_generalized .+= torque_third_law
    return nothing
end


# function addGeneralizedForcesThirdLaw!(wrench::Wrench{T}, tm::TypedMechanismScenario{N,T}, cInfo::ContactInstructions) where {N,T}
#     addGeneralizedForcesExternal!( wrench, tm, tm.bodyBodyCache.mesh_1.BodyID)
#     addGeneralizedForcesExternal!(-wrench, tm, tm.bodyBodyCache.mesh_2.BodyID)
#     return nothing
# end
#
# function addGeneralizedForcesExternal!(wrench::Wrench{T}, tm::TypedMechanismScenario{N,T}, body_id::BodyID) where {N,T}
#     jac = tm.GeometricJacobian[body_id]
#     if jac != nothing
#         torque!(tm.torque_third_law, jac, wrench)
#         tm.f_generalized .+= tm.torque_third_law
#     end
#     return nothing
# end



# function calcXd!(xx::AbstractVector{T1}, x::AbstractVector{T1}, m::MechanismScenario{NX,NQ,T2},
#         tm::TypedMechanismScenario{NQ,T1}, t::Float64=0.0) where {NX,NQ,T1,T2}
#
#     state = tm.state
#     copyto!(tm, x)
#     H = tm.result.massmatrix
#     mass_matrix!(H, state)
#     dynamics_bias!(tm.result, state)
#     configuration_derivative!(tm.result.q̇, state)
#     forceAllElasticIntersections!(m, tm)
#
#     (m.continuous_controller == nothing) || m.continuous_controller(tm, t)
#
#     f_generalized = tm.f_generalized
#     rhs = tm.result.dynamicsbias.parent
#     rhs .*= -1.0
#     rhs .+= f_generalized
#     rhs .+= m.τ_ext
#     rhs .+= tm.τ_ext
#
#     chol_fact = LinearAlgebra.cholesky!(H)
#     ldiv!(tm.result.v̇.parent, chol_fact, rhs)
#
#     copyto!(xx, tm, tm.result)
#     return nothing
# end



# f_generalized = tm.f_generalized
# rhs = tm.result.dynamicsbias.parent
# rhs .*= -1.0
# rhs .+= f_generalized
# rhs .+= m.τ_ext
# rhs .+= tm.τ_ext
#
# chol_fact = LinearAlgebra.cholesky!(H)
# ldiv!(tm.result.v̇.parent, chol_fact, rhs)



#
