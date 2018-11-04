
function verify_bristle_ids!(m::MechanismScenario{NX,NQ,T}, x::Vector{Float64}) where {NX,NQ,T}
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

function calcXd!(xx::AbstractVector{T}, x::AbstractVector{T}, m::MechanismScenario{NX,NQ,T}) where {NX,NQ,T}
    return calcXd!(xx, x, m, m.dual)
end

function calcXd!(xx::AbstractVector{Float64}, x::AbstractVector{Float64}, m::MechanismScenario{NX,NQ,T}) where {NX,NQ,T}
    return calcXd!(xx, x, m, m.float)
end

function calcXd!(xx::AbstractVector{T1}, x::AbstractVector{T1}, m::MechanismScenario{NX,NQ,T2},
        tm::TypedMechanismScenario{NQ,T1}) where {NX,NQ,T1,T2}

    state = tm.state
    # any(isnan.(x)) && error("nan in x")
    copyto!(tm, x)
    # any(isnan.(state.q)) && error("nan in state.q")  # isnan
    # any(isnan.(state.v)) && error("nan in state.v")  # isnan
    H = tm.result.massmatrix
    mass_matrix!(H, state)
    dynamics_bias!(tm.result, state)
    configuration_derivative!(tm.result.q̇, state)
    forceAllElasticIntersections!(m, tm)
    f_generalized = tm.f_generalized
    # any(isnan.(f_generalized)) && error("nan in tm.f_gen_cum")  # isnan
    # C = tm.result.dynamicsbias.parent
    # rhs = -C + f_generalized + m.τ_ext
    rhs = tm.result.dynamicsbias.parent
    rhs .*= -1.0
    rhs .+= f_generalized
    rhs .+= m.τ_ext
    #
    # tm.result.v̇ .= H \ rhs
    chol_fact = LinearAlgebra.cholesky!(H)
    ldiv!(tm.result.v̇.parent, chol_fact, rhs)
    #
    copyto!(xx, tm, tm.result)
    return nothing
end

function forceAllElasticIntersections!(m::MechanismScenario{NX,NQ,T1}, tm::TypedMechanismScenario{NQ,T2}) where {NX,NQ,T1,T2}
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
                integrate_over_logic!(tm.bodyBodyCache, m.TT_Cache)
                if !isempty(tm.bodyBodyCache.TractionCache)
                    is_no_wrench = false
                    if is_bristle
                        wrench = bristle_friction!(m.frame_world, tm, con_ins_k)
                    else
                        wrench = regularized_friction(m.frame_world, tm.bodyBodyCache)
                    end
                    tm.bodyBodyCache.wrench = wrench
                    addGeneralizedForcesThirdLaw!(wrench, tm, con_ins_k)
                end
            end
        end
        (is_bristle & is_no_wrench) && bristle_friction_no_contact!(tm, con_ins_k)
    end
    return nothing
end

function refreshJacobians!(m::MechanismScenario{NX,NQ,T1}, tm::TypedMechanismScenario{NQ,T2}) where {NX,NQ,T1,T2}
    for k = m.body_ids
        path_k = m.path[k]
        (path_k != nothing) && geometric_jacobian!(tm.GeometricJacobian[k], tm.state, path_k)
    end
    return nothing
end

function normal_wrench(frame::CartesianFrame3D, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    wrench = zero(Wrench{T}, frame)
    @inbounds begin
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            p_dA = calc_point_p_dA(trac, k)
            wrench += Wrench(trac.r_cart[k], p_dA * trac.n̂)
        end
    end
    end
    return wrench
end

function calcTriTetIntersections!(m::MechanismScenario, con_ins_k::ContactInstructions) # where {N,T}
    b = m.float.bodyBodyCache  # this can be float because intersection is assumed to not depend on partials
    refreshBodyBodyTransform!(m, m.float, con_ins_k)
    x_r¹_r² = inv(b.x_rʷ_r¹) * b.x_rʷ_r²
    update_TT_Cache!(m.TT_Cache, translation(x_r¹_r²), rotation(x_r¹_r²))
    if con_ins_k.mutual_compliance
        tree_tree_intersect(m.TT_Cache, b.mesh_1.tet.tet.tree, b.mesh_2.tet.tet.tree)
    else
        tree_tree_intersect(m.TT_Cache, b.mesh_1.tri.tree, b.mesh_2.tet.tet.tree)
    end
    return nothing
end

function refreshBodyBodyTransform!(m::MechanismScenario, tm::TypedMechanismScenario{N,T},
        con_ins_k::ContactInstructions) where {N,T}

    b = tm.bodyBodyCache
    b.mesh_1 = m.MeshCache[con_ins_k.id_1]
    b.mesh_2 = m.MeshCache[con_ins_k.id_2]
    b.x_rʷ_r¹ = transform_to_root(tm.state, b.mesh_1.BodyID)  # TODO: add safe=false
    b.x_rʷ_r² = transform_to_root(tm.state, b.mesh_2.BodyID)  # TODO: add safe=false
    b.x_r²_rʷ = inv(b.x_rʷ_r²)
    b.x_r¹_rʷ = inv(b.x_rʷ_r¹)
    return nothing
end

function refreshBodyBodyCache!(m::MechanismScenario, tm::TypedMechanismScenario{N,T},
        con_ins_k::ContactInstructions) where {N,T}

    b = tm.bodyBodyCache
    empty!(b.TractionCache)
    refreshBodyBodyTransform!(m, tm, con_ins_k)

    twist_w_r¹ = twist_wrt_world(tm.state, b.mesh_1.BodyID)
    twist_w_r² = twist_wrt_world(tm.state, b.mesh_2.BodyID)
    b.twist_r¹_r² = -twist_w_r² + twist_w_r¹  # velocity of tri wrt tet exp in world

    b.μ = con_ins_k.μ_pair
    mat_2 = b.mesh_2.tet.c_prop
    b.Ē = mat_2.Ē
    b.χ = mat_2.χ
    b.d⁻¹ = mat_2.d⁻¹
    return nothing
end

function integrate_over_logic!(b::TypedElasticBodyBodyCache{N,T}, ttCache::TT_Cache) where {N,T}
    mesh_1 = b.mesh_1
    mesh_2 = b.mesh_2
    x_rʷ_r¹ = b.x_rʷ_r¹
    x_rʷ_r² = b.x_rʷ_r²
    x_r¹_rʷ = b.x_r¹_rʷ
    x_r²_rʷ = b.x_r²_rʷ
    if is_compliant(mesh_1)
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
    ind_vert = m.tri.ind[i_tri]
    return m.point[ind_vert]
end

function tetrahedron_vertices_ϵ(i_tet::Int64, m::MeshCache)
    ind_vert = m.tet.tet.ind[i_tet]
    ϵ = m.tet.ϵ[ind_vert]
    ϵ = SMatrix{1,4,Float64,4}(ϵ)
    cart_vert = m.point[ind_vert]
    return cart_vert, ϵ
end

function calc_ζ_transforms(frame_ζ::CartesianFrame3D, frame_r ::CartesianFrame3D, p_tet, x_r_w, x_w_r)
    x_r_ζ = MatrixTransform(frame_ζ, frame_r, asMatOnePad(p_tet))
    x_w_ζ = x_w_r * x_r_ζ
    x_ζ_w = inv(x_r_ζ) * x_r_w  # NOTE: inv(A_r¹_ζ) is **always** Float64
    return x_w_ζ, x_ζ_w
end

function find_plane_tet(E::Float64, ϵ::SMatrix{1,4,Float64,4}, X_r_w)
    return (E * ϵ) * X_r_w
end

function integrate_over_volume_volume!(i_1::Int64, i_2::Int64, mesh_1::MeshCache, mesh_2::MeshCache,
        x_rʷ_r¹::Transform3D{T}, x_rʷ_r²::Transform3D{T}, x_r¹_rʷ::Transform3D{T}, x_r²_rʷ::Transform3D{T},
        b::TypedElasticBodyBodyCache{N,T}) where {N,T}

    vert_1, ϵ¹ = tetrahedron_vertices_ϵ(i_1, mesh_1)
    vert_2, ϵ² = tetrahedron_vertices_ϵ(i_2, mesh_2)
    x_rʷ_ζ¹, x_ζ¹_rʷ = calc_ζ_transforms(FRAME_ζ¹, mesh_1.FrameID, vert_1, x_r¹_rʷ, x_rʷ_r¹)
    x_rʷ_ζ², x_ζ²_rʷ = calc_ζ_transforms(FRAME_ζ², mesh_2.FrameID, vert_2, x_r²_rʷ, x_rʷ_r²)

    plane_1_rʷ = find_plane_tet(get_Ē(mesh_1), ϵ¹, x_ζ¹_rʷ.mat)
    plane_2_rʷ = find_plane_tet(get_Ē(mesh_2), ϵ², x_ζ²_rʷ.mat)
    plane_rʷ = plane_2_rʷ - plane_1_rʷ  # pressure_2 - pressure_1

    poly_rʷ = clip_plane_tet(plane_rʷ, x_rʷ_ζ¹.mat)
    if 3 <= length(poly_rʷ)
        poly_ζ² = one_pad_then_mul(x_ζ²_rʷ.mat, poly_rʷ)
        poly_ζ² = zero_small_coordinates(poly_ζ²)  # This needs to be done to avoid a degenerate situation where the
            # plane lies exactly on the intersection of the faces of two tets. This situation happens EVERY time two
            # tet faces that lie on the surface intersect.
        poly_ζ² = clip_in_tet_coordinates(poly_ζ²)
        if 3 <= length(poly_ζ²)
            frame_world = b.frame_world
            n = unPad(plane_rʷ)
            n̂ = unsafe_normalize(n)
            n̂ = FreeVector3D(frame_world, n̂)
            integrate_over_polygon_patch!(b, poly_ζ², frame_world, n̂, x_rʷ_ζ², x_ζ²_rʷ, ϵ²)
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
    x_rʷ_ζ², x_ζ²_rʷ = calc_ζ_transforms(FRAME_ζ², mesh_2.FrameID, vert_2, x_r²_rʷ, x_rʷ_r²)
    x_ζ²_r¹ = x_ζ²_rʷ * x_rʷ_r¹

    n̂_r¹ = FreeVector3D(mesh_1.FrameID, -triangleNormal(vert_1))  # pressure is applied opposite the trianle normal
    n̂_rʷ = x_rʷ_r¹ * n̂_r¹

    poly_rʷ = poly_eight(vert_1)
    poly_ζ² = one_pad_then_mul(x_ζ²_r¹.mat, poly_rʷ)
    poly_ζ² = clip_in_tet_coordinates(poly_ζ²)
    if 3 <= length(poly_ζ²)
        frame_world = b.frame_world
        integrate_over_polygon_patch!(b, poly_ζ², frame_world, n̂_rʷ, x_rʷ_ζ², x_ζ²_rʷ, ϵ²)
    end
end

function integrate_over_polygon_patch!(b::TypedElasticBodyBodyCache{N,T}, poly_ζ²::poly_eight{4,T},
        frame_world::CartesianFrame3D, n̂::FreeVector3D{SVector{3,T}}, x_rʷ_ζ²::MatrixTransform{4,4,T,16},
        x_ζ²_rʷ::MatrixTransform{4,4,T,16}, ϵ::SMatrix{1,4,Float64,4}) where {N,T}

    poly_rʷ = mul_then_un_pad(x_rʷ_ζ².mat, poly_ζ²)
    centroid_rʷ = Point3D(frame_world, centroid(poly_rʷ)[2])
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
        area_quad_k = area(vert_rʷ_1, vert_rʷ_2, centroid_rʷ)
        (0.0 < area_quad_k) && fillTractionCacheForTriangle!(b, area_quad_k, x_ζ²_ϕ, x_rʷ_ζ², n̂, ϵ)  # no point in adding intersection if area is zero
    end
    return nothing
end

# # TODO: create fillTractionCacheForTriangle! with a macro
function fillTractionCacheForTriangle!(b::TypedElasticBodyBodyCache{1,T}, area_quad_k::T,
        A_ζ_ϕ::MatrixTransform{4,3,T,12}, A_w_ζ::MatrixTransform{4,4,T,16}, n̂::FreeVector3D{SVector{3,T}},
        ϵ::SMatrix{1,4,Float64,4}) where {T}

    r_cart_1, v_cart_t_1, pene_1, dA_1, p_1 = fillTractionCacheInnerLoop!(1, b, A_ζ_ϕ, A_w_ζ, n̂, area_quad_k, ϵ)
    p = (p_1, )
    if sum(p) != 0.0
        dA = (dA_1, )
        r_cart = (r_cart_1, )
        v_cart_t = (v_cart_t_1, )
        pene = (pene_1, )
        trac_cache = TractionCache(n̂, r_cart, v_cart_t, pene, dA, p)
        addCacheItem!(b.TractionCache, trac_cache)
    end
    return nothing
end

function fillTractionCacheForTriangle!(b::TypedElasticBodyBodyCache{3,T}, area_quad_k::T,
        A_ζ_ϕ::MatrixTransform{4,3,T,12}, A_w_ζ::MatrixTransform{4,4,T,16}, n̂::FreeVector3D{SVector{3,T}},
        ϵ::SMatrix{1,4,Float64,4}) where {T}

    r_cart_1, v_cart_t_1, pene_1, dA_1, p_1 = fillTractionCacheInnerLoop!(1, b, A_ζ_ϕ, A_w_ζ, n̂, area_quad_k, ϵ)
    r_cart_2, v_cart_t_2, pene_2, dA_2, p_2 = fillTractionCacheInnerLoop!(2, b, A_ζ_ϕ, A_w_ζ, n̂, area_quad_k, ϵ)
    r_cart_3, v_cart_t_3, pene_3, dA_3, p_3 = fillTractionCacheInnerLoop!(3, b, A_ζ_ϕ, A_w_ζ, n̂, area_quad_k, ϵ)
    p = (p_1, p_2, p_3)
    if sum(p) != 0.0
        dA = (dA_1, dA_2, dA_3)
        r_cart = (r_cart_1, r_cart_2, r_cart_3)
        v_cart_t = (v_cart_t_1, v_cart_t_2, v_cart_t_3)
        pene = (pene_1, pene_2, pene_3)
        trac_cache = TractionCache(n̂, r_cart, v_cart_t, pene, dA, p)
        addCacheItem!(b.TractionCache, trac_cache)
    end
    return nothing
end

function fillTractionCacheInnerLoop!(k::Int64, b::TypedElasticBodyBodyCache{N,T},
        A_ζ_ϕ::MatrixTransform{4,3,T,12}, A_w_ζ::MatrixTransform{4,4,T,16}, n̂::FreeVector3D{SVector{3,T}},
        area_quad_k::T, ϵ::SMatrix{1,4,Float64,4}) where {N,T}

    quad_point_ϕ = Point3D(FRAME_ϕ, b.quad.zeta[k])
    quad_point_ζ = A_ζ_ϕ * quad_point_ϕ
    p_cart_qp = unPad(A_w_ζ * quad_point_ζ)
    ϵ_quad = dot(ϵ, quad_point_ζ.v)
    cart_vel, signed_mag_vel_n = calcNormalVelocityMag(b.twist_r¹_r², p_cart_qp, n̂)
    ϵ_dot = signed_mag_vel_n * b.d⁻¹  # ϵ_dot ≈ z_dot / thickness because the rigid body provides the normal and is fixed
    damp_term = fastSoftPlus(1.0 - b.χ * ϵ_dot)
    p = -ϵ_quad * b.Ē
    p_hc = p * damp_term  # -ϵ_quad * damp_term * b.Ē
    dA = b.quad.w[k] * area_quad_k
    penetration = ϵ_quad / b.d⁻¹   # normal penetration (l = p L / E)
    return p_cart_qp, cart_vel, penetration, dA, p_hc
end

function calcNormalVelocityMag(twist_tri_tet::Twist{T}, p_cart_qp::Point3D{SVector{3,T}}, n̂::FreeVector3D{SVector{3,T}}) where {T}
    cart_vel = point_velocity(twist_tri_tet, p_cart_qp)
    signed_mag_vel_n = dot(n̂, cart_vel)
    return cart_vel, signed_mag_vel_n
end

function calcTangentialVelocity(twist_tri_tet::Twist{T}, p_cart_qp::Point3D{SVector{3,T}}, n̂::FreeVector3D{SVector{3,T}}) where {T}
    cart_vel_crw = point_velocity(twist_tri_tet, p_cart_qp)
    signed_mag_vel_n = dot(n̂, cart_vel_crw)
    cart_vel_crw_n = n̂ * signed_mag_vel_n
    cart_vel_crw_t = cart_vel_crw - cart_vel_crw_n
    return cart_vel_crw_t, signed_mag_vel_n
end

function addGeneralizedForcesThirdLaw!(wrench::Wrench{T}, tm::TypedMechanismScenario{N,T}, cInfo::ContactInstructions) where {N,T}
    addGeneralizedForcesExternal!(wrench,  1.0, tm, tm.bodyBodyCache.mesh_1.BodyID)
    addGeneralizedForcesExternal!(wrench, -1.0, tm, tm.bodyBodyCache.mesh_2.BodyID)
    return nothing
end

function addGeneralizedForcesExternal!(wrench::Wrench{T},  coeff::Float64, tm::TypedMechanismScenario{N,T},
        body_id::BodyID) where {N,T}

    jac = tm.GeometricJacobian[body_id]
    if jac != nothing
        tm.f_generalized .+= coeff * torque(jac, wrench)
    end
    return nothing
end
