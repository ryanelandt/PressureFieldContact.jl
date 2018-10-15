# function regularized_friction(frame::CartesianFrame3D, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
#     wrench = zeroWrench(frame, T)
#     for k_trac = 1:length(b.TractionCache)
#         trac = b.TractionCache[k_trac]
#         for k = 1:N
#             cart_vel_crw_t = trac.v_cart_t[k]
#             mag_vel_t  = norm(cart_vel_crw_t.v)
#             μ_reg     = b.μ * fastSigmoid(mag_vel_t)
#             traction_k = trac.p_dA[k] * (trac.traction_normal - μ_reg * safe_normalize(cart_vel_crw_t))
#             wrench += Wrench(trac.r_cart[k], traction_k)
#         end
#     end
#     return wrench
# end

function regularized_friction(frame::CartesianFrame3D, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    wrench = zero(Wrench{T}, frame)

    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            cart_vel_crw_t = trac.v_cart_t[k]
            mag_vel_t = safe_norm(cart_vel_crw_t.v)
            μ_reg = b.μ * fastSigmoid(mag_vel_t)
            traction_k = trac.p_dA[k] * (trac.traction_normal - μ_reg * safe_normalize(cart_vel_crw_t))
            wrench += Wrench(trac.r_cart[k], traction_k)
        end
    end
    return wrench
end

function find_contact_pressure_center(b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    int_p_dA = zero(T)
    int_p_r_dA = zeros(SVector{3,T})
    @inbounds begin
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            p = trac.p_dA[k]
            int_p_dA += p
            int_p_r_dA += p * trac.r_cart[k].v
        end
    end
    end
    frame = b.TractionCache[1].r_cart[1].frame
    return Point3D(frame, int_p_r_dA / int_p_dA ), int_p_dA
end

function normal_wrench(frame::CartesianFrame3D, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    wrench = zeroWrench(frame, T)
    @inbounds begin
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            wrench += Wrench(trac.r_cart[k], trac.p_dA[k] * trac.traction_normal)
        end
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

function calc_T_θ_dA(b::TypedElasticBodyBodyCache{N,T}, τ_θ_s_w::FreeVector3D{SVector{3,T}}, r̄_w::Point3D{SVector{3,T}}) where {N,T}
    @framecheck(τ_θ_s_w.frame, r̄_w.frame)
    int_T_θ_dA = zero(T)
    @inbounds begin
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            r_rel_w = trac.r_cart[k] - r̄_w
            r_cross_τ = cross(r_rel_w, τ_θ_s_w)
            r_cross_τ_squared = norm_squared(r_cross_τ)
            int_T_θ_dA += trac.p_dA[k] * r_cross_τ_squared
        end
    end
    end
    return int_T_θ_dA
end

function notch_function(x::T, k=0.05) where {T}
    fact_k = pi / k
    if abs(fact_k * x) < Float64(pi)
        return 0.5 * (1 - cos(fact_k * x))
    else
        return one(T)
    end
end

function bristle_friction_inner(b::TypedElasticBodyBodyCache{N,T}, BF::BristleFriction, c_ins::ContactInstructions,
    frame_w::CartesianFrame3D, c_θ::FreeVector3D{SVector{3,T}}, c_r::FreeVector3D{SVector{3,T}}) where {N,T}

    λ_r_s_w, τ_θ_s_w, r̄_w, p_dA_patch = bristle_no_slip_force_moment(b, frame_w, BF, c_θ, c_r)
    T_θ_dA = calc_T_θ_dA(b, τ_θ_s_w, r̄_w)
    wrench_λ_r̄_w = zeroWrench(frame_w, T)
    const_λ_term = λ_r_s_w * safe_scalar_divide(one(T), p_dA_patch)
    @inbounds begin
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        n̂_w = trac.traction_normal
        for k = 1:N
            p = trac.p_dA[k]
            r_rel_w = trac.r_cart[k] - r̄_w
            r_cross_τ = cross(r_rel_w, τ_θ_s_w)
            term_θ_num = norm_squared(r_cross_τ)
            term_θ_den = norm_squared(r_rel_w) * T_θ_dA
            λ_goal_w = const_λ_term - r_cross_τ * safe_scalar_divide(term_θ_num, term_θ_den)
            σ_f = vec_sub_vec_proj(λ_goal_w, n̂_w)

            notch_ratio = safe_scalar_divide(norm_squared(σ_f), norm_squared(λ_goal_w))
            notch_factor = notch_function(notch_ratio)

            σ̂_f = safe_normalize(σ_f) * notch_factor
            traction_t_w = σ̂_f * p * c_ins.μ_pair
            wrench_λ_r̄_w += Wrench(cross(r_rel_w, traction_t_w), traction_t_w)
        end
    end
    end

    λ_r = FreeVector3D(frame_w, linear(wrench_λ_r̄_w))
    τ_θ = FreeVector3D(frame_w, angular(wrench_λ_r̄_w))
    λ_r = soft_clamp(λ_r, λ_r_s_w)
    τ_θ = soft_clamp(τ_θ, τ_θ_s_w)
    ang_term = τ_θ + cross(r̄_w, λ_r)
    wrench_λ_w = Wrench(frame_w, ang_term.v, λ_r.v)
    λ_r = transform(λ_r, b.x_tet_root)
    τ_θ = transform(τ_θ, b.x_tet_root)

    return wrench_λ_w, τ_θ, λ_r
end


function bristle_no_slip_force_moment(b::TypedElasticBodyBodyCache{N,T}, frame_w::CartesianFrame3D,
    BF::BristleFriction, c_θ::FreeVector3D{SVector{3,T}}, c_r::FreeVector3D{SVector{3,T}}) where {N,T}

    inv_τ = 1 / BF.τ
    twist_r_c_w = b.twist_tri_tet
    r̄_w, p_dA_patch = find_contact_pressure_center(b)
    ċ_r_rel = point_velocity(twist_r_c_w, r̄_w)
    c_w_rel = FreeVector3D(frame_w, angular(twist_r_c_w))
    c_θ = transform(c_θ, b.x_root_tet)
    c_r = transform(c_r, b.x_root_tet)
    λ_r_s_w = -BF.K_r * (c_r + inv_τ * ċ_r_rel)
    τ_θ_s_w = -BF.K_θ * (c_θ + inv_τ * c_w_rel)
    return λ_r_s_w, τ_θ_s_w, r̄_w, p_dA_patch
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
    wrench_λ_w, τ_θ, λ_r = bristle_friction_inner(b, BF, c_ins, frame_w, c_θ, c_r)
    wrench_t_w = normal_wrench(frame_w, b) + wrench_λ_w
    c_w = -BF.τ * (c_θ + τ_θ * (1 / BF.K_θ) )
    ċ_r = -BF.τ * (c_r + λ_r * (1 / BF.K_r) )
    update_bristle_d1!(tm, bristle_id, c_w, ċ_r)
    return wrench_t_w
end
