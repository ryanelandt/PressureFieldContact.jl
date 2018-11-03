
function regularized_friction(frame::CartesianFrame3D, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    wrench = zero(Wrench{T}, frame)

    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            cart_vel = trac.v_cart[k]
            n̂ = trac.n̂
            cart_vel_t = vec_sub_vec_proj(cart_vel, n̂)
            mag_vel_t = safe_norm(cart_vel_t.v)
            μ_reg = b.μ * fastSigmoid(mag_vel_t)
            p_dA = calc_point_p_dA(trac, k)
            traction_k = p_dA * (n̂ - μ_reg * safe_normalize(cart_vel_t))
            wrench += Wrench(trac.r_cart[k], traction_k)
        end
    end
    return wrench
end

function bristle_friction_no_contact!(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    BF = c_ins.BristleFriction
    bristle_id = BF.BristleID
    segments_s = get_bristle_d0(tm, bristle_id)
    segments_ṡ = get_bristle_d1(tm, bristle_id)
    segments_ṡ .= -BF.τ * segments_s
    return nothing
end

function bristle_friction!(frame_w::CartesianFrame3D, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    b = tm.bodyBodyCache
    BF = c_ins.BristleFriction
    bristle_id = BF.BristleID
    frame_c = b.mesh_2.FrameID
    r̄, int_p_dA = find_contact_pressure_center(b::TypedElasticBodyBodyCache{N,T})
    K_θ = bristle_stiffness_θ_now(b, r̄, frame_c)
    K_r = bristle_stiffness_l_now(b, frame_c)
    c_θ, c_r = bristle_deformation_warp(frame_c, tm, bristle_id, K_θ, K_r)
    wrench_λ_w, τ_θ, λ_r = bristle_friction_inner(b, BF, c_ins, frame_w, c_θ, c_r, K_θ, K_r)
    wrench_t_w = normal_wrench(frame_w, b) + wrench_λ_w
    c_w_e, ċ_r_e = bristle_deformation_rate(tm, bristle_id, K_θ, K_r, τ_θ, λ_r, BF, frame_c)
    update_bristle_d1_warp!(tm, bristle_id, c_w_e, ċ_r_e)
    return wrench_t_w
end

function bristle_stiffness_l_now(b::TypedElasticBodyBodyCache{N,T}, frame_c::CartesianFrame3D) where {N,T}
    K_l_dimless = zeros(SVector{3,T})
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        n̂ = trac.n̂.v
        for k = 1:N
            dA = trac.dA[k]
            K_l_dA_1 = safe_norm(vec_sub_vec_proj(SVector{3,T}(1.0, 0.0, 0.0), n̂))
            K_l_dA_2 = safe_norm(vec_sub_vec_proj(SVector{3,T}(0.0, 1.0, 0.0), n̂))
            K_l_dA_3 = safe_norm(vec_sub_vec_proj(SVector{3,T}(0.0, 0.0, 1.0), n̂))
            K_l_dimless += dA * SVector{3,T}(K_l_dA_1, K_l_dA_2, K_l_dA_3)
        end
    end
    K_r = b.Ē * K_l_dimless * b.mesh_2.tet.c_prop.d⁻¹
    all(0.0 .<= K_r) || error("element negative")
    K_r = FreeVector3D(b.frame_world, K_r)
    return transform(K_r, b.x_r²_rʷ)
end

function bristle_stiffness_θ_now(b::TypedElasticBodyBodyCache{N,T}, r̄::Point3D{SVector{3,T}}, frame_c::CartesianFrame3D) where {N,T}

    K_θ_dimless = zeros(SVector{3,T})
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        n̂ = trac.n̂.v
        for k = 1:N
            dA = trac.dA[k]
            r_rel = trac.r_cart[k] - r̄
            norm_r_rel = dot(r_rel.v, r_rel.v)  # safe_norm(r_rel.v)
            K_θ_dA_1 = abs(dot(n̂, SVector{3,Float64}(1.0, 0.0, 0.0)))
            K_θ_dA_2 = abs(dot(n̂, SVector{3,Float64}(0.0, 1.0, 0.0)))
            K_θ_dA_3 = abs(dot(n̂, SVector{3,Float64}(0.0, 0.0, 1.0)))
            K_θ_dimless += dA * norm_r_rel * SVector{3,T}(K_θ_dA_1, K_θ_dA_2, K_θ_dA_3)
        end
    end
    K_θ = b.Ē * K_θ_dimless * b.mesh_2.tet.c_prop.d⁻¹
    all(0.0 .<= K_θ) || error("element negative")
    K_θ = FreeVector3D(b.frame_world, K_θ)
    return transform(K_θ, b.x_r²_rʷ)
end

function bristle_deformation_e(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID) where {N,T}
    segments_s = get_bristle_d0(tm, bristle_id)
    c_θ_ref = SVector{3,T}(segments_s[1], segments_s[2], segments_s[3])
    c_l_ref = translation(SPQuatFloating{Float64}(), segments_s)
    return c_θ_ref, c_l_ref
end

function bristle_deformation_warp(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID, K_θ::SVector{3,T},
        K_r::SVector{3,T}) where {N,T}

    c_θ_ref, c_r_ref = bristle_deformation_e(tm, bristle_id)
    c_θ = c_θ_ref ./ sqrt.(K_θ)
    c_r = c_r_ref ./ sqrt.(K_r)
    return c_θ, c_r
end

function bristle_deformation_warp(frame_c::CartesianFrame3D, tm::TypedMechanismScenario{N,T}, bristle_id::BristleID,
        K_θ::FreeVector3D{SVector{3,T}}, K_r::FreeVector3D{SVector{3,T}}) where {N,T}

    @framecheck(K_r.frame, frame_c)
    c_θ, c_r = bristle_deformation_warp(tm, bristle_id, K_r.v, K_θ.v)
    c_r = FreeVector3D(frame_c, c_r)
    c_θ = FreeVector3D(frame_c, c_θ)
    return c_θ, c_r
end

function bristle_deformation_rate(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID,
        K_θ::FreeVector3D{SVector{3,T}}, K_r::FreeVector3D{SVector{3,T}}, τ_θ::FreeVector3D{SVector{3,T}},
        λ_r::FreeVector3D{SVector{3,T}}, BF::BristleFriction, frame_c::CartesianFrame3D) where {N,T}
    #

    c_θ, c_r = bristle_deformation_warp(frame_c, tm, bristle_id, K_θ, K_r)

    @framecheck(c_θ.frame, τ_θ.frame)
    @framecheck(c_θ.frame, K_θ.frame)
    @framecheck(c_θ.frame, frame_c)
    c_w = -BF.τ * (c_θ.v + τ_θ.v ./ K_θ.v)

    @framecheck(c_r.frame, λ_r.frame)
    @framecheck(c_r.frame, K_r.frame)
    @framecheck(c_r.frame, frame_c)
    ċ_r = -BF.τ * (c_r.v + λ_r.v ./ K_r.v)

    c_θ_e, c_r_e = bristle_deformation_e(tm, bristle_id)
    c_θθ_e = sqrt.(K_θ.v) .* c_w
    c_rr_e = sqrt.(K_r.v) .* ċ_r
    return c_θθ_e, c_rr_e
end

calc_point_p_dA(tc::TractionCache{N,T}, k::Int64) where {N,T} = tc.p[k] * tc.dA[k]

function find_contact_pressure_center(b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    int_p_dA = zero(T)
    int_p_r_dA = zeros(SVector{3,T})
    @inbounds begin
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            p = calc_point_p_dA(trac, k)
            int_p_dA += p
            int_p_r_dA += p * trac.r_cart[k].v
        end
    end
    end
    frame = b.TractionCache[1].r_cart[1].frame
    return Point3D(frame, int_p_r_dA / int_p_dA ), int_p_dA
end

@inline get_bristle_d0(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID) where {N,T} = segments(tm.s)[bristle_id]
@inline get_bristle_d1(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID) where {N,T} = segments(tm.ṡ)[bristle_id]

function update_bristle_d1_warp!(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID, c_w::FreeVector3D{SVector{3,T}},
        ċ_r::FreeVector3D{SVector{3,T}}) where {N,T}

    @framecheck(c_w.frame, ċ_r.frame)
    update_bristle_d1_warp!(tm, bristle_id, c_w.v, ċ_r.v)
    return nothing
end

function update_bristle_d1_warp!(tm::TypedMechanismScenario{N,T}, bristle_id::BristleID, c_w::SVector{3,T},
        ċ_r::SVector{3,T}) where {N,T}

    segments_ṡ = get_bristle_d1(tm, bristle_id)
    segments_ṡ .= vcat(c_w, ċ_r)
    return nothing
end

function calc_T_θ_dA(b::TypedElasticBodyBodyCache{N,T}, τ_θ_s_w::FreeVector3D{SVector{3,T}},
        r̄_w::Point3D{SVector{3,T}}) where {N,T}

    @framecheck(τ_θ_s_w.frame, r̄_w.frame)
    int_T_θ_dA = zero(T)
    @inbounds begin
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            r_rel_w = trac.r_cart[k] - r̄_w
            r_cross_τ = cross(r_rel_w, τ_θ_s_w)
            r_cross_τ_squared = norm_squared(r_cross_τ)
            p_dA = calc_point_p_dA(trac, k)
            int_T_θ_dA += p_dA * r_cross_τ_squared
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
        frame_w::CartesianFrame3D, c_θ::FreeVector3D{SVector{3,T}}, c_r::FreeVector3D{SVector{3,T}},
        K_θ::FreeVector3D{SVector{3,T}}, K_r::FreeVector3D{SVector{3,T}}) where {N,T}

    λ_r_s_w, τ_θ_s_w, r̄_w, p_dA_patch = bristle_no_slip_force_moment(b, frame_w, BF, c_θ, c_r, K_θ, K_r)
    T_θ_dA = calc_T_θ_dA(b, τ_θ_s_w, r̄_w)
    wrench_λ_r̄_w = zero(Wrench{T}, frame_w)
    const_λ_term = λ_r_s_w * safe_scalar_divide(one(T), p_dA_patch)
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        n̂_w = trac.n̂
        for k = 1:N
            r_rel_w = trac.r_cart[k] - r̄_w
            r_cross_τ = cross(r_rel_w, τ_θ_s_w)
            term_θ_num = norm_squared(r_cross_τ)
            term_θ_den = norm_squared(r_rel_w) * T_θ_dA
            λ_goal_w = const_λ_term - r_cross_τ * safe_scalar_divide(term_θ_num, term_θ_den)
            σ_f = vec_sub_vec_proj(λ_goal_w, n̂_w)
            notch_ratio = safe_scalar_divide(norm_squared(σ_f), norm_squared(λ_goal_w))
            notch_factor = notch_function(notch_ratio, 0.15)
            σ̂_f = safe_normalize(σ_f) * notch_factor
            p_dA = calc_point_p_dA(trac, k)
            traction_t_w = σ̂_f * p_dA * c_ins.μ_pair
            wrench_λ_r̄_w += Wrench(cross(r_rel_w, traction_t_w), traction_t_w)
        end
    end

    λ_r = FreeVector3D(frame_w, linear(wrench_λ_r̄_w))
    τ_θ = FreeVector3D(frame_w, angular(wrench_λ_r̄_w))
    λ_r = soft_clamp(λ_r, λ_r_s_w)
    τ_θ = soft_clamp(τ_θ, τ_θ_s_w)
    ang_term = τ_θ + cross(r̄_w, λ_r)
    wrench_λ_w = Wrench(frame_w, ang_term.v, λ_r.v)
    λ_r = transform(λ_r, b.x_r²_rʷ)
    τ_θ = transform(τ_θ, b.x_r²_rʷ)

    return wrench_λ_w, τ_θ, λ_r
end

function bristle_no_slip_force_moment(b::TypedElasticBodyBodyCache{N,T}, frame_w::CartesianFrame3D,
        BF::BristleFriction, c_θ::FreeVector3D{SVector{3,T}}, c_r::FreeVector3D{SVector{3,T}},
        K_θ::FreeVector3D{SVector{3,T}}, K_r::FreeVector3D{SVector{3,T}}) where {N,T}

    inv_τ = 1 / BF.τ
    twist_r¹_r²_rʷ = b.twist_r¹_r²

    r̄_w, p_dA_patch = find_contact_pressure_center(b)
    ċ_r_rel = point_velocity(twist_r¹_r²_rʷ, r̄_w)
    c_w_rel = FreeVector3D(frame_w, angular(twist_r¹_r²_rʷ))
    c_θ = transform(c_θ, b.x_rʷ_r²)
    c_r = transform(c_r, b.x_rʷ_r²)
    K_θ = transform(K_θ, b.x_rʷ_r²)
    K_r = transform(K_r, b.x_rʷ_r²)
    @framecheck(K_θ.frame, c_θ.frame)
    τ_θ_s_w = (c_θ + inv_τ * c_w_rel)
    τ_θ_s_w = FreeVector3D(K_θ.frame, -K_θ.v .* τ_θ_s_w.v)
    @framecheck(K_r.frame, c_r.frame)
    λ_r_s_w = (c_r + inv_τ * ċ_r_rel)
    λ_r_s_w = FreeVector3D(K_r.frame, -K_r.v .* λ_r_s_w.v)
    return λ_r_s_w, τ_θ_s_w, r̄_w, p_dA_patch
end
