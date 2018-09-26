calcXd!(xx::AbstractVector{T}, x::AbstractVector{T}, m::MechanismScenario{N,T}) where {N,T} = calcXd!(xx, x, m, m.dual)
calcXd!(xx::AbstractVector{Float64}, x::AbstractVector{Float64}, m::MechanismScenario{N,T}) where {N,T} = calcXd!(xx, x, m, m.float)
function calcXd!(xx::AbstractVector{T1}, x::AbstractVector{T1}, m::MechanismScenario{N,T2}, tm::TypedMechanismScenario{N,T1}) where {N,T1,T2}
  state = tm.state
  copyto!(state, x)
  # any(isnan.(state.q)) && error("nan in state.q")  # isnan
  # any(isnan.(state.v)) && error("nan in state.v")  # isnan
  H = tm.result.massmatrix
  mass_matrix!(H, state)
  dynamics_bias!(tm.result, state)
  configuration_derivative!(tm.result.q̇, state)
  forceAllElasticIntersections!(m, tm)
  f_generalized = tm.f_generalized
  # any(isnan.(f_generalized)) && error("nan in tm.f_gen_cum")  # isnan
  C = tm.result.dynamicsbias.parent
  rhs = -C + f_generalized + m.tau_ext
  tm.result.v̇ .= H \ rhs
  copyto!(xx, tm.result)
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
        if 0 != length(m.TT_Cache)
            refreshBodyBodyCache!(m, tm, con_ins_k)
            integrateOverContactPatch!(tm.bodyBodyCache, m.TT_Cache)
            wrench_zero = Wrench(Point3D(m.frame_world, SVector{3,T2}(0.0, 0.0, 0.0)), FreeVector3D(m.frame_world, SVector{3,T2}(0.0, 0.0, 0.0)))
            wrench = friction_model(con_ins_k.friction_model, wrench_zero, tm.bodyBodyCache)
            addGeneralizedForcesThirdLaw!(wrench, tm, con_ins_k)
        end
    end
    return nothing
end

function calcTriTetIntersections!(m::MechanismScenario, con_ins_k::ContactInstructions) where {N,T}
    b = m.float.bodyBodyCache  # this can be float because intersection is assumed to not depend on partials
    refreshBodyBodyTransform!(m, m.float, con_ins_k)
    x_tri_tet = inv(b.x_root_tri) * b.x_root_tet
    update_TT_Cache!(m.TT_Cache, translation(x_tri_tet), rotation(x_tri_tet))
    tree_tree_intersect(m.TT_Cache, b.mesh_tri.raw.tri_tree, b.mesh_tet.raw.tet_tree)
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

    mat_tet = b.mesh_tet.raw.material
    twist_root_tri = twist_wrt_world(tm.state, b.mesh_tri.BodyID)
    twist_root_tet = twist_wrt_world(tm.state, b.mesh_tet.BodyID)

    b.twist_tri_tet = -twist_root_tet + twist_root_tri  # velocity of tri wrt tet exp in world

    b.frac_epsilon = con_ins_k.frac_epsilon
    b.mu = con_ins_k.mu_pair

    b.E_effective = mat_tet.contact.E_effective
    b.hc_velocity_damping = mat_tet.contact.hc_velocity_damping
    b.inv_thickness = mat_tet.contact.inv_thickness
    return nothing
end

function integrateOverContactPatch!(b::TypedElasticBodyBodyCache{N,T}, ttCache::TT_Cache) where {N,T}
    mesh_tri = b.mesh_tri
    mesh_tet = b.mesh_tet
    triTetCache = b.triTetCache
    quadTriCache = triTetCache.quadTriCache
    for k = 1:length(ttCache.vc)
        i_tri, i_tet = ttCache.vc[k]
        p_tri = mesh_tri.raw.point[mesh_tri.raw.tri_ind[i_tri]]
        i_vert_tet = mesh_tet.raw.tet_ind[i_tet]
        p_tet = mesh_tet.raw.point[i_vert_tet]
        triTetCache.strain = mesh_tet.raw.strain[i_vert_tet]
        # tet calculation
        A_tet_zeta = MatrixTransform(frame_tet_cs, mesh_tet.FrameID, asMatOnePad(p_tet))
        triTetCache.A_w_zeta_top = getTop(b.x_root_tet * A_tet_zeta)
        A_zeta_w = inv(A_tet_zeta) * b.x_tet_root  # NOTE: inv(A_tet_zeta) is **always** Float64
        # tri calculation
        triTetCache.traction_normal = b.x_root_tri * FreeVector3D(mesh_tri.FrameID, -triangleNormal(p_tri))
        A_zeta_r = A_zeta_w * b.x_root_tri
        poly_4D_1 = quadTriCache.clip_poly_4D_1
        poly_4D_2 = quadTriCache.clip_poly_4D_2
        poly_3D = quadTriCache.clip_poly_3D
        empty!(poly_4D_1)
        for j = 1:3
            add!(poly_4D_1, A_zeta_r * Point3D(mesh_tri.FrameID, p_tri[j]))
        end
        modularTriTetClip(poly_4D_1, poly_4D_2)
        if 3 <= poly_4D_1.i  ### Clip ###
            tet_clip_poly_to_cartesian!(poly_3D, poly_4D_1, triTetCache.A_w_zeta_top)
            poly_area, centroid_w_no_frame = centroid(poly_3D)
            triTetCache.centroid_w = Point3D(poly_3D.frame, centroid_w_no_frame)  # TODO: add to Tri_Tet_Intersections
            if 0.0 < poly_area
                triTetCache.centroid_zeta = A_zeta_w * triTetCache.centroid_w
                intersectionTriangulation!(b)
            end
        end
    end
    return nothing
end

function intersectionTriangulation!(b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    triTetCache = b.triTetCache
    quadTriCache = triTetCache.quadTriCache
    poly_4D = quadTriCache.clip_poly_4D_1
    poly_3D = quadTriCache.clip_poly_3D

    n_vertices = poly_4D.i
    zeta_2 = getPoint(poly_4D, n_vertices)
    vert_tri_2 = getPoint(poly_3D, n_vertices)
    for k_cent_tri = 1:n_vertices
        zeta_1 = zeta_2
        zeta_2 = getPoint(poly_4D, k_cent_tri)
        quadTriCache.A_zeta_phi = hcat(zeta_1, zeta_2, triTetCache.centroid_zeta)

        vert_tri_1 = vert_tri_2
        vert_tri_2 = getPoint(poly_3D, k_cent_tri)
        quadTriCache.area_quad_k = area(vert_tri_1, vert_tri_2, triTetCache.centroid_w)
        (0.0 < quadTriCache.area_quad_k) && fillTractionCacheForTriangle!(b)  # no point in adding intersection if area is zero
    end
    return nothing
end

function fillTractionCacheForTriangle!(b::TypedElasticBodyBodyCache{N,T}) where {N,T}
  trac_now = returnNext(b.TractionCache)
  trac_now.traction_normal = b.triTetCache.traction_normal
  for k = 1:N
    fillTractionCacheInnerLoop!(k, trac_now, b)
  end
  all(trac_now.p_dA .== 0.0) && (b.TractionCache.ind_fill -= 1)  # pressure was zero at all points
  return nothing
end

function fillTractionCacheInnerLoop!(k::Int64, trac_now::TractionCache{N,T}, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    triTetCache = b.triTetCache
    quadTriCache = triTetCache.quadTriCache

    quad_point_phi = Point3D(frame_tri_cs, b.quad.zeta[k])
    quad_point_zeta = quadTriCache.A_zeta_phi * quad_point_phi
    p_cart_qp = triTetCache.A_w_zeta_top * quad_point_zeta
    strain_quad = b.frac_epsilon * dot(triTetCache.strain, quad_point_zeta.v)
    cart_vel_crw_t, signed_mag_vel_n = calcTangentialVelocity(b.twist_tri_tet, p_cart_qp, triTetCache.traction_normal)
    epsilon_dot = b.frac_epsilon * signed_mag_vel_n * b.inv_thickness  # epsilon_dot ≈ z_dot / thickness because the rigid body provides the normal and is fixed
    damp_term = fastSoftPlus(1.0 - b.hc_velocity_damping * epsilon_dot)
    the_pressure = -strain_quad * damp_term

    trac_now.r_cart[k]   = p_cart_qp
    trac_now.v_cart_t[k] = cart_vel_crw_t
    trac_now.p_dA[k]     = the_pressure * b.quad.w[k] * quadTriCache.area_quad_k * b.E_effective
    return nothing
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
