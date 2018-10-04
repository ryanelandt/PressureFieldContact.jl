function regularized_friction(frame::CartesianFrame3D, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    wrench = zeroWrench(frame, T)
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            cart_vel_crw_t = trac.v_cart_t[k]
            mag_vel_t  = norm(cart_vel_crw_t.v)
            mu_reg     = b.mu * fastSigmoid(mag_vel_t)
            traction_k = trac.p_dA[k] * (trac.traction_normal - mu_reg * safe_normalize(cart_vel_crw_t))
            wrench += Wrench(trac.r_cart[k], traction_k)
        end
    end
    return wrench
end

function find_contact_pressure_center(b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    int_p_dA = zero(T)
    int_p_r_dA = zeros(SVector{3,T})
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            p = trac.p_dA[k]
            int_p_dA += p
            int_p_r_dA += p * trac.r_cart[k].v
        end
    end
    frame = b.TractionCache[1].r_cart[1].frame
    return Point3D(frame, int_p_r_dA / int_p_dA ), int_p_dA
end

function normal_wrench(frame::CartesianFrame3D, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    wrench = zeroWrench(frame, T)
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            wrench += Wrench(trac.r_cart[k], trac.p_dA[k] * trac.traction_normal)
        end
    end
    return wrench
end

function bristle_deformation(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID) where {N,T}
    segments_s = get_bristle_d0(tm, bristle_id)
    c_θ = rotation(SPQuatFloating{Float64}(), segments_s)
    c_θ = cheapRV(c_θ)
    c_r = translation(SPQuatFloating{Float64}(), segments_s)
    return c_θ, c_r
end

function bristle_deformation(frame_c::CartesianFrame3D, tm::TypedMechanismScenario{N,T}, bristle_id::BristleID) where {N,T}
    c_θ, c_r = bristle_deformation(tm, bristle_id)
    c_r = FreeVector3D(frame_c, c_r)
    c_θ = FreeVector3D(frame_c, c_θ)
    return c_θ, c_r
end

@inline get_bristle_d0(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID) where {N,T} = segments(tm.s)[bristle_id]
@inline get_bristle_d1(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID) where {N,T} = segments(tm.ṡ)[bristle_id]
function update_bristle_d1!(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID, c_w::FreeVector3D{SVector{3,T}}, ċ_r::FreeVector3D{SVector{3,T}}) where {N,T}
    @framecheck(c_w.frame, ċ_r.frame)
    update_bristle_d1!(tm, bristle_id, c_w.v, ċ_r.v)
    return nothing
end

function update_bristle_d1!(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID, c_w::SVector{3,T}, ċ_r::SVector{3,T}) where {N,T}
    segments_s = get_bristle_d0(tm, bristle_id)
    segments_ṡ = get_bristle_d1(tm, bristle_id)
    velocity_to_configuration_derivative!(segments_ṡ, SPQuatFloating{Float64}(), segments_s, vcat(c_w, ċ_r))
    return nothing
end

function calc_λ_θ_s(r_rel_c::FreeVector3D{SVector{3,T}}, τ_θ_s::FreeVector3D{SVector{3,T}}) where {T}
    return cross(r_rel_c, τ_θ_s) * (-safe_inv_norm_squared(r_rel_c.v))
end

function calc_T_θ(r_rel_c::FreeVector3D{SVector{3,T}}, τ_θ_s::FreeVector3D{SVector{3,T}}, p::T) where {T}
    v2 = cross(r_rel_c, τ_θ_s)
    return p * dot(v2, v2)
end

function calc_T_θ_dA(b::TypedElasticBodyBodyCache{N,T}, τ_θ_s::FreeVector3D, r̄_c::Point3D) where {N,T}
    @framecheck(τ_θ_s.frame, r̄_c.frame)
    int_T_θ_dA = zero(T)
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            r_cart_c = transform(trac.r_cart[k], b.x_tet_root)
            r_rel_c = r_cart_c - r̄_c
            int_T_θ_dA += calc_T_θ(r_rel_c, τ_θ_s, trac.p_dA[k])
        end
    end
    return int_T_θ_dA
end

function bristle_friction_inner(b::TypedElasticBodyBodyCache{N,T}, BF, c_ins, frame_c, c_θ, c_r) where {N,T}
    λ_r_s, τ_θ_s, r̄_c, p_dA_patch = bristle_no_slip_force_moment(b, frame_c, BF, c_θ, c_r)
    wrench_λ_c = zeroWrench(frame_c, T)
    T_θ_dA = calc_T_θ_dA(b, τ_θ_s, r̄_c)
    τ_θ_s_w = transform(τ_θ_s, b.x_root_tet)
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        n̂_c = transform(trac.traction_normal, b.x_tet_root)
        for k = 1:N
            p = trac.p_dA[k]
            r_c = transform(trac.r_cart[k], b.x_tet_root)
            r_rel_c = r_c - r̄_c
            W_r = safe_scalar_divide(p, p_dA_patch)
            W_θ = safe_scalar_divide(calc_T_θ(r_rel_c, τ_θ_s, p), T_θ_dA)
            λ_goal_c = W_r * λ_r_s + W_θ * calc_λ_θ_s(r_rel_c, τ_θ_s)
            σ_f = vec_sub_vec_proj(λ_goal_c, n̂_c)
            σ̂_f = safe_normalize(σ_f)
            traction_t_c = σ̂_f * p * c_ins.mu_pair
            wrench_λ_c += Wrench(r_c, traction_t_c)
        end
    end
    λ_r = FreeVector3D(frame_c, linear(wrench_λ_c))
    τ_θ = FreeVector3D(frame_c, angular(wrench_λ_c)) - cross(r̄_c, λ_r)
    λ_r = soft_clamp(λ_r, λ_r_s)
    τ_θ = soft_clamp(τ_θ, τ_θ_s)
    ang_term = τ_θ + cross(r̄_c, λ_r)
    wrench_λ_c = Wrench(frame_c, ang_term.v, λ_r.v)
    wrench_λ_w = transform(wrench_λ_c, b.x_root_tet)
    return wrench_λ_w, τ_θ, λ_r
end

function bristle_no_slip_force_moment(b::TypedElasticBodyBodyCache{N,T}, frame_c, BF::BristleFriction, c_θ, c_r) where {N,T}
    inv_τ = 1 / BF.τ
    twist_r_c_w = b.twist_tri_tet
    twist_r_c_c = transform(twist_r_c_w, b.x_tet_root)
    r̄_w, p_dA_patch = find_contact_pressure_center(b)
    r̄_c = transform(r̄_w, b.x_tet_root)
    ċ_r_rel = point_velocity(twist_r_c_c, r̄_c)
    c_w_rel = FreeVector3D(frame_c, angular(twist_r_c_c))
    λ_r_s = -BF.K_r * (c_r + inv_τ * ċ_r_rel)
    τ_θ_s = -BF.K_θ * (c_θ + inv_τ * c_w_rel);
    return λ_r_s, τ_θ_s, r̄_c, p_dA_patch
end

function bristle_friction_no_contact!(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    BF = c_ins.BristleFriction
    bristle_id = BF.BristleID
    c_θ, c_r = bristle_deformation(tm, bristle_id)
    c_w = -BF.τ * c_θ
    ċ_r = -BF.τ * c_r
    update_bristle_d1!(tm, bristle_id, c_w, ċ_r)
    return nothing
end

function bristle_friction!(frame_w::CartesianFrame3D, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    b = tm.bodyBodyCache
    BF = c_ins.BristleFriction
    bristle_id = BF.BristleID
    frame_c = b.mesh_tet.FrameID
    c_θ, c_r = bristle_deformation(frame_c, tm, bristle_id)
    wrench_λ_w, τ_θ, λ_r = bristle_friction_inner(b, BF, c_ins, frame_c, c_θ, c_r)
    wrench_n̂_w = normal_wrench(frame_w, b)
    wrench_t_w = wrench_n̂_w + wrench_λ_w
    c_w = -BF.τ * (c_θ + τ_θ * (1 / BF.K_θ) )
    ċ_r = -BF.τ * (c_r + λ_r * (1 / BF.K_r) )
    update_bristle_d1!(tm, bristle_id, c_w, ċ_r)
    return wrench_t_w
end
