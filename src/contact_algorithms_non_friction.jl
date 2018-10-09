function verify_bristle_ids!(m::MechanismScenario{N,T}, x::Vector{Float64}) where {N,T}
    copyto!(m.float, x)

    bristle_done = falses(m.bristle_ids)
    for k = 1:length(m.ContactInstructions)
        con_ins_k = m.ContactInstructions[k]
        if con_ins_k.BristleFriction != nothing
            bristle_id = con_ins_k.BristleFriction.BristleID
            (bristle_done[bristle_id] == true) && error("BristleID $bristle_id assigned twice")
            bristle_done[bristle_id] = true
        end
    end
    for bristle_id = m.bristle_ids
        if !bristle_done[bristle_id]
            segments_s = get_bristle_d0(m.float, bristle_id)
            fill!(segments_s, 0.0)
            println("BristleID $bristle_id not assigned")
        end
    end

    copyto!(x, m.float)
    return nothing
end

calcXd!(xx::AbstractVector{T}, x::AbstractVector{T}, m::MechanismScenario{N,T}) where {N,T} = calcXd!(xx, x, m, m.dual)
calcXd!(xx::AbstractVector{Float64}, x::AbstractVector{Float64}, m::MechanismScenario{N,T}) where {N,T} = calcXd!(xx, x, m, m.float)
function calcXd!(xx::AbstractVector{T1}, x::AbstractVector{T1}, m::MechanismScenario{N,T2}, tm::TypedMechanismScenario{N,T1}) where {N,T1,T2}
    state = tm.state
    any(isnan.(x)) && error("nan in x")
    copyto!(tm, x)
    any(isnan.(state.q)) && error("nan in state.q")  # isnan
    any(isnan.(state.v)) && error("nan in state.v")  # isnan
    H = tm.result.massmatrix
    mass_matrix!(H, state)
    dynamics_bias!(tm.result, state)
    configuration_derivative!(tm.result.q̇, state)
    forceAllElasticIntersections!(m, tm)
    f_generalized = tm.f_generalized
    any(isnan.(f_generalized)) && error("nan in tm.f_gen_cum")  # isnan
    C = tm.result.dynamicsbias.parent
    rhs = -C + f_generalized + m.τ_ext
    tm.result.v̇ .= H \ rhs
    copyto!(xx, tm, tm.result)
    return nothing
end

function refreshJacobians!(m::MechanismScenario{N,T1}, tm::TypedMechanismScenario{N,T2}) where {N,T1,T2}
    for k = m.body_ids
        path_k = m.path[k]
        (path_k != nothing) && geometric_jacobian!(tm.GeometricJacobian[k], tm.state, path_k)
    end
    return nothing
end

function forceAllElasticIntersections!(m::MechanismScenario{N,T1}, tm::TypedMechanismScenario{N,T2}) where {N,T1,T2}
    refreshJacobians!(m, tm)
    tm.f_generalized .= zero(T2)
    for k = 1:length(m.ContactInstructions)
        con_ins_k = m.ContactInstructions[k]
        calcTriTetIntersections!(m, con_ins_k)
        is_intersections = (0 != length(m.TT_Cache))
        is_bristle = (con_ins_k.BristleFriction != nothing)
        is_no_wrench = true
        if is_intersections | is_bristle
            if is_intersections
                refreshBodyBodyCache!(m, tm, con_ins_k)
                integrateOverContactPatch!(tm.bodyBodyCache, m.TT_Cache)
                if !isempty(tm.bodyBodyCache.TractionCache)
                    is_no_wrench = false
                    if is_bristle
                        wrench = bristle_friction!(m.frame_world, tm, con_ins_k)
                    else
                        wrench = regularized_friction(m.frame_world, tm.bodyBodyCache)
                    end
                    addGeneralizedForcesThirdLaw!(wrench, tm, con_ins_k)
                end
            end
        end
        (is_bristle & is_no_wrench) && bristle_friction_no_contact!(tm, con_ins_k)
    end
    return nothing
end

function calcTriTetIntersections!(m::MechanismScenario, con_ins_k::ContactInstructions) where {N,T}
    b = m.float.bodyBodyCache  # this can be float because intersection is assumed to not depend on partials
    refreshBodyBodyTransform!(m, m.float, con_ins_k)
    x_tri_tet = inv(b.x_root_tri) * b.x_root_tet
    update_TT_Cache!(m.TT_Cache, translation(x_tri_tet), rotation(x_tri_tet))
    tree_tree_intersect(m.TT_Cache, b.mesh_tri.tri.tree, b.mesh_tet.tet.tet.tree)
    return nothing
end

function refreshBodyBodyTransform!(m::MechanismScenario, tm::TypedMechanismScenario{N,T}, con_ins_k::ContactInstructions) where {N,T}
    b = tm.bodyBodyCache
    b.mesh_tri = m.MeshCache[con_ins_k.id_tri]
    b.mesh_tet = m.MeshCache[con_ins_k.id_tet]
    b.x_root_tri = transform_to_root(tm.state, b.mesh_tri.BodyID)  # TODO: add safe=false
    b.x_root_tet = transform_to_root(tm.state, b.mesh_tet.BodyID)  # TODO: add safe=false
    b.x_tet_root = inv(b.x_root_tet)
    return nothing
end

function refreshBodyBodyCache!(m::MechanismScenario, tm::TypedMechanismScenario{N,T}, con_ins_k::ContactInstructions) where {N,T}
    b = tm.bodyBodyCache
    empty!(b.TractionCache)
    refreshBodyBodyTransform!(m, tm, con_ins_k)

    mat_tet = b.mesh_tet.tet.contact_prop
    twist_root_tri = twist_wrt_world(tm.state, b.mesh_tri.BodyID)
    twist_root_tet = twist_wrt_world(tm.state, b.mesh_tet.BodyID)

    b.twist_tri_tet = -twist_root_tet + twist_root_tri  # velocity of tri wrt tet exp in world

    b.frac_ϵ = con_ins_k.frac_ϵ
    b.μ = con_ins_k.μ_pair

    b.Ē = mat_tet.Ē
    b.χ = mat_tet.χ
    b.d⁻¹ = mat_tet.d⁻¹
    return nothing
end

function integrateOverContactPatch!(b::TypedElasticBodyBodyCache{N,T}, ttCache::TT_Cache) where {N,T}
    mesh_tri = b.mesh_tri
    mesh_tet = b.mesh_tet
    triTetCache = b.triTetCache
    quadTriCache = triTetCache.quadTriCache
    x_root_tet = b.x_root_tet
    x_tet_root = b.x_tet_root
    x_root_tri = b.x_root_tri
    for k = 1:length(ttCache.vc)
        i_tri, i_tet = ttCache.vc[k]
        p_tri = mesh_tri.point[mesh_tri.tri.ind[i_tri]]
        i_vert_tet = mesh_tet.tet.tet.ind[i_tet]
        p_tet = mesh_tet.point[i_vert_tet]
        triTetCache.ϵ = mesh_tet.tet.ϵ[i_vert_tet]
        # tet calculation
        A_tet_ζ = MatrixTransform(frame_tet_cs, mesh_tet.FrameID, asMatOnePad(p_tet))
        # triTetCache.A_w_ζ_top = getTop(x_root_tet * A_tet_ζ)
        A_w_ζ_top = getTop(x_root_tet * A_tet_ζ)
        A_ζ_w = inv(A_tet_ζ) * x_tet_root  # NOTE: inv(A_tet_ζ) is **always** Float64
        # tri calculation
        triTetCache.traction_normal = x_root_tri * FreeVector3D(mesh_tri.FrameID, -triangleNormal(p_tri))
        A_ζ_r = A_ζ_w * x_root_tri
        poly_4D_1 = quadTriCache.clip_poly_4D_1
        poly_4D_2 = quadTriCache.clip_poly_4D_2
        poly_3D = quadTriCache.clip_poly_3D
        empty!(poly_4D_1)
        for j = 1:3
            add!(poly_4D_1, A_ζ_r * Point3D(mesh_tri.FrameID, p_tri[j]))
        end
        modularTriTetClip(poly_4D_1, poly_4D_2)
        if 3 <= poly_4D_1.i  ### Clip ###
            # tet_clip_poly_to_cartesian!(poly_3D, poly_4D_1, triTetCache.A_w_ζ_top)
            tet_clip_poly_to_cartesian!(poly_3D, poly_4D_1, A_w_ζ_top)
            poly_area, centroid_w_no_frame = centroid(poly_3D)
            triTetCache.centroid_w = Point3D(poly_3D.frame, centroid_w_no_frame)  # TODO: add to Tri_Tet_Intersections
            if 0.0 < poly_area
                triTetCache.centroid_ζ = A_ζ_w * triTetCache.centroid_w
                intersectionTriangulation!(b, A_w_ζ_top)
            end
        end
    end
    return nothing
end

function intersectionTriangulation!(b::TypedElasticBodyBodyCache{N,T}, A_w_ζ_top) where {N,T}
    triTetCache = b.triTetCache
    quadTriCache = triTetCache.quadTriCache
    poly_4D = quadTriCache.clip_poly_4D_1
    poly_3D = quadTriCache.clip_poly_3D

    n_vertices = poly_4D.i
    ζ_2 = getPoint(poly_4D, n_vertices)
    vert_tri_2 = getPoint(poly_3D, n_vertices)
    for k_cent_tri = 1:n_vertices
        ζ_1 = ζ_2
        ζ_2 = getPoint(poly_4D, k_cent_tri)
        # quadTriCache.A_ζ_ϕ = hcat(ζ_1, ζ_2, triTetCache.centroid_ζ)
        A_ζ_ϕ = hcat(ζ_1, ζ_2, triTetCache.centroid_ζ)

        vert_tri_1 = vert_tri_2
        vert_tri_2 = getPoint(poly_3D, k_cent_tri)
        quadTriCache.area_quad_k = area(vert_tri_1, vert_tri_2, triTetCache.centroid_w)
        (0.0 < quadTriCache.area_quad_k) && fillTractionCacheForTriangle!(b, A_ζ_ϕ, A_w_ζ_top)  # no point in adding intersection if area is zero
    end
    return nothing
end

# function fillTractionCacheForTriangle!(b::TypedElasticBodyBodyCache{N,T}) where {N,T}
#   trac_now = returnNext(b.TractionCache)
#   trac_now.traction_normal = b.triTetCache.traction_normal
#   for k = 1:N
#     fillTractionCacheInnerLoop!(k, trac_now, b)
#   end
#   all(trac_now.p_dA .== 0.0) && (b.TractionCache.ind_fill -= 1)  # pressure was zero at all points
#   return nothing
# end
#
# function fillTractionCacheInnerLoop!(k::Int64, trac_now::TractionCache{N,T}, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
#     triTetCache = b.triTetCache
#     quadTriCache = triTetCache.quadTriCache
#
#     quad_point_ϕ = Point3D(frame_tri_cs, b.quad.zeta[k])
#     quad_point_ζ = quadTriCache.A_ζ_ϕ * quad_point_ϕ
#     p_cart_qp = triTetCache.A_w_ζ_top * quad_point_ζ
#     ϵ_quad = b.frac_ϵ * dot(triTetCache.ϵ, quad_point_ζ.v)
#     cart_vel_crw_t, signed_mag_vel_n = calcTangentialVelocity(b.twist_tri_tet, p_cart_qp, triTetCache.traction_normal)
#     ϵ_dot = b.frac_ϵ * signed_mag_vel_n * b.d⁻¹  # ϵ_dot ≈ z_dot / thickness because the rigid body provides the normal and is fixed
#     damp_term = fastSoftPlus(1.0 - b.χ * ϵ_dot)
#     p = -ϵ_quad * damp_term
#
#     trac_now.r_cart[k]   = p_cart_qp
#     trac_now.v_cart_t[k] = cart_vel_crw_t
#     trac_now.p_dA[k]     = p * b.quad.w[k] * quadTriCache.area_quad_k * b.Ē
#     return nothing
# end

# struct TractionCache{N,T}
#     traction_normal::FreeVector3D{SVector{3,T}}
#     r_cart::NTuple{N,Point3D{SVector{3,T}}}
#     v_cart_t::NTuple{N,FreeVector3D{SVector{3,T}}}
#     p_dA::NTuple{N,T}
#     function TractionCache{N,T}(frame::CartesianFrame3D) where {N,T}
#         traction_normal = FreeVector3D(frame, SVector{3, T}(NaN .+ zeros(3)))
#         r_cart = NTuple{N,Point3D{SVector{3,Float64}}}([Point3D(frame, SVector{3,Float64}(NaN, NaN, NaN)) for k = 1:3])
#         v_cart_t = NTuple{N,FreeVector3D{SVector{3,Float64}}}([FreeVector3D(frame, SVector{3,Float64}(NaN, NaN, NaN)) for k = 1:3])
#         p_dA = NTuple{N, T}(NaN .+ zeros(N))
#         return new(traction_normal, r_cart, v_cart_t, p_dA)
#     end
# end


# struct TractionCache{N,T}
#     traction_normal::FreeVector3D{SVector{3,T}}
#     r_cart::NTuple{N,Point3D{SVector{3,T}}}
#     v_cart_t::NTuple{N,FreeVector3D{SVector{3,T}}}
#     p_dA::NTuple{N,T}
#     function TractionCache{N,T}(traction_normal::FreeVector3D{SVector{3,T}}, r_cart::NTuple{N,Point3D{SVector{3,T}}},
#                                 v_cart_t::NTuple{N,FreeVector3D{SVector{3,T}}}, p_dA::NTuple{N,T})
#
#        return new(traction_normal, r_cart, v_cart_t, p_dA)
#     end
# end


function fillTractionCacheForTriangle!(b::TypedElasticBodyBodyCache{3,T}, A_ζ_ϕ, A_w_ζ_top) where {T}
  # trac_now = returnNext(b.TractionCache)
  # trac_now.traction_normal =
  traction_normal = b.triTetCache.traction_normal
  r_cart_1, v_cart_t_1, p_dA_1 = fillTractionCacheInnerLoop!(1, b, A_ζ_ϕ, A_w_ζ_top)
  r_cart_2, v_cart_t_2, p_dA_2 = fillTractionCacheInnerLoop!(2, b, A_ζ_ϕ, A_w_ζ_top)
  r_cart_3, v_cart_t_3, p_dA_3 = fillTractionCacheInnerLoop!(3, b, A_ζ_ϕ, A_w_ζ_top)
  p_dA = (p_dA_1, p_dA_2, p_dA_3)
  if sum(p_dA) != 0.0
      r_cart = (r_cart_1, r_cart_2, r_cart_3)
      v_cart_t = (v_cart_t_1, v_cart_t_2, v_cart_t_3)
      trac_cache = TractionCache(traction_normal, r_cart, v_cart_t, p_dA)
      addCacheItem!(b.TractionCache, trac_cache)
  end

  # for k = 1:N
    # fillTractionCacheInnerLoop!(k, trac_now, b)
  # end
  # all(trac_now.p_dA .== 0.0) && (b.TractionCache.ind_fill -= 1)  # pressure was zero at all points
  return nothing
end

function fillTractionCacheInnerLoop!(k::Int64, b::TypedElasticBodyBodyCache{N,T}, A_ζ_ϕ, A_w_ζ_top) where {N,T}
    triTetCache = b.triTetCache
    quadTriCache = triTetCache.quadTriCache

    quad_point_ϕ = Point3D(frame_tri_cs, b.quad.zeta[k])
    # quad_point_ζ = quadTriCache.A_ζ_ϕ * quad_point_ϕ
    # p_cart_qp = triTetCache.A_w_ζ_top * quad_point_ζ
    quad_point_ζ = A_ζ_ϕ * quad_point_ϕ
    p_cart_qp = A_w_ζ_top * quad_point_ζ
    #
    ϵ_quad = b.frac_ϵ * dot(triTetCache.ϵ, quad_point_ζ.v)
    cart_vel_crw_t, signed_mag_vel_n = calcTangentialVelocity(b.twist_tri_tet, p_cart_qp, triTetCache.traction_normal)
    ϵ_dot = b.frac_ϵ * signed_mag_vel_n * b.d⁻¹  # ϵ_dot ≈ z_dot / thickness because the rigid body provides the normal and is fixed
    damp_term = fastSoftPlus(1.0 - b.χ * ϵ_dot)
    p = -ϵ_quad * damp_term
    # trac_now.r_cart[k]   = p_cart_qp
    # trac_now.v_cart_t[k] = cart_vel_crw_t
    # trac_now.p_dA[k]     = p * b.quad.w[k] * quadTriCache.area_quad_k * b.Ē
    p_dA = p * b.quad.w[k] * quadTriCache.area_quad_k * b.Ē
    return p_cart_qp, cart_vel_crw_t, p_dA
end

function calcTangentialVelocity(twist_tri_tet::Twist{T}, p_cart_qp::Point3D{SVector{3,T}}, traction_normal::FreeVector3D{SVector{3,T}}) where {T}
    cart_vel_crw = point_velocity(twist_tri_tet, p_cart_qp)
    signed_mag_vel_n = dot(traction_normal, cart_vel_crw)
    cart_vel_crw_n = traction_normal * signed_mag_vel_n
    cart_vel_crw_t = cart_vel_crw - cart_vel_crw_n
    return cart_vel_crw_t, signed_mag_vel_n
end

function addGeneralizedForcesThirdLaw!(wrench::Wrench{T}, tm::TypedMechanismScenario{N,T}, cInfo::ContactInstructions) where {N,T}
    coeff = cInfo.frac_linear_weight
    addGeneralizedForcesExternal!(wrench,  coeff, tm, tm.bodyBodyCache.mesh_tri.BodyID)
    addGeneralizedForcesExternal!(wrench, -coeff, tm, tm.bodyBodyCache.mesh_tet.BodyID)
    return nothing
end

function addGeneralizedForcesExternal!(wrench::Wrench{T},  coeff::Float64, tm::TypedMechanismScenario{N,T}, body_id::BodyID) where {N,T}
    jac = tm.GeometricJacobian[body_id]
    if jac != nothing
        tm.f_generalized .+= coeff * torque(jac, wrench)
    end
    return nothing
end
