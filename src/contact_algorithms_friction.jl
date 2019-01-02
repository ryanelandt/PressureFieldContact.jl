
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
            p_dA = calc_p_dA(trac, k)
            traction_k = p_dA * (n̂ - μ_reg * safe_normalize(cart_vel_t))
            wrench += Wrench(trac.r_cart[k], traction_k)
        end
    end
    return wrench
end

function normal_wrench_patch_center(frame::CartesianFrame3D, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    wrench = zero(Wrench{T}, frame)
    int_p_r_dA = zeros(SVector{3,T})
    int_p_dA = zero(T)
    # @inbounds begin
        for k_trac = 1:length(b.TractionCache)
            trac = b.TractionCache[k_trac]
            for k = 1:N
                p_dA = calc_p_dA(trac, k)
                wrench += Wrench(trac.r_cart[k], p_dA * trac.n̂)
                int_p_r_dA += trac.r_cart[k].v * p_dA
                int_p_dA += p_dA
            end
        end
    # end
    p_center = Point3D(frame, int_p_r_dA / int_p_dA)
    return wrench, p_center
end

function make_stiffness_PD!(K::MMatrix{6,6,T,36}) where {T}
    # NOTE: # diag(SMatrix{6,6,Int64,36}(collect(1:36)))

    tol = 1.0e-4
    ang_K = tol * (K[1] + K[8] + K[15])
    K[1] += ang_K
    K[8] += ang_K
    K[15] += ang_K
    lin_K = tol * (K[22] + K[29] + K[36])
    K[22] += lin_K
    K[29] += lin_K
    K[36] += lin_K
end

function veil_friction_no_contact!(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    BF = c_ins.BristleFriction
    bristle_id = BF.BristleID
    segments_s = get_bristle_d0(tm, bristle_id)
    segments_ṡ = get_bristle_d1(tm, bristle_id)
    segments_ṡ .= -BF.τ * segments_s
    return nothing
end

##########################

function calc_patch_spatial_stiffness(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions, K) where {N,T}

    # TODO: do this in a frame rotated with the world frame and coincident with the body frame

    b = tm.bodyBodyCache
    frameᶜ = b.mesh_2.FrameID
    tc = b.TractionCache

    Q_11_sum = zeros(SMatrix{3,3,T,9})
    Q_12_sum = zeros(SMatrix{3,3,T,9})
    Q_22_sum = zeros(SMatrix{3,3,T,9})
    for k = 1:length(tc)
        trac = tc.vec[k]
        n̂² = transform(trac.n̂, b.x_r²_rʷ)
        Q_22 = I - n̂².v * n̂².v'  # suprisingly fast
        for k_qp = 1:N
            r² = transform(trac.r_cart[k_qp], b.x_r²_rʷ)
            r²_skew = vector_to_skew_symmetric(r².v)
            Q_12 = r²_skew * Q_22
            minus_Q_11 = Q_12 * r²_skew
            p_dA = calc_p_dA(trac, k_qp)
            Q_11_sum -= p_dA * minus_Q_11
            Q_12_sum += p_dA * Q_12
            Q_22_sum += p_dA * Q_22
        end
    end
    μ = c_ins.μ_pair
    K[1:3, 1:3] .= μ * Q_11_sum
    K[1:3, 4:6] .= μ * Q_12_sum
    K[4:6, 4:6] .= μ * Q_22_sum
    return nothing
end

function veil_friction!(frameʷ::CartesianFrame3D, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    # TODO: do all calculations relative to patch center

    BF = c_ins.BristleFriction
    bristle_id = BF.BristleID
    Δ = get_bristle_d0(tm, bristle_id)
    Δ = SVector{6,T}(Δ[1], Δ[2], Δ[3], Δ[4], Δ[5], Δ[6])

    b = tm.bodyBodyCache
    twist_r¹_r² = b.twist_r¹_r²
    twist_r¹_r²_r² = transform(twist_r¹_r², b.x_r²_rʷ)

    frameᶜ = b.mesh_2.FrameID
    @framecheck(twist_r¹_r²_r².frame, frameᶜ)
    v_spatial_rel = as_static_vector(twist_r¹_r²_r²)

    K = b.K
    fill!(K.data, zero(T))
    calc_patch_spatial_stiffness(tm, c_ins, K.data)

    make_stiffness_PD!(K.data)
    K.data .*= BF.k̄

    K_err = Matrix(K)
    if isposdef(K_err)
        #
    else
        println("")
        for k_33 = 1:6
            println(K_err[k_33, 1:6])
        end
        println("")
        println("smoking_gun = [")
        for k = 1:36
            a4 = K[k]
            @printf("%.17f", a4)
            (k == 36) || (print(","))
        end
        println("]")
        println("")
        println("length(b.TractionCache): ", length(b.TractionCache))
        println("")
        error("something is very very wrong")
    end

    # NOTE:
    # Forward       Inverse
    # K = Uᵀ U      K⁻¹ = U⁻¹ U⁻ᵀ
    # Δ = U δ       δ = U⁻¹ Δ
    # ΔΔ = U δδ     δδ = U⁻¹ ΔΔ

    cf = cholesky(K)
    U = cf.U
    U⁻¹ = SMatrix{6,6,T,36}(inv(U))
    δ = U⁻¹ * Δ
    wrench²_un = calc_spatial_bristle_force(tm, c_ins, frameᶜ, twist_r¹_r²_r², δ)

    ######################

    τ_s2 = -K * (Δ + BF.τ * v_spatial_rel)
    wrench_stick_2 = Wrench{T}(frameᶜ, SVector{3,T}(τ_s2[1], τ_s2[2], τ_s2[3]), SVector{3,T}(τ_s2[4], τ_s2[5], τ_s2[6]))
    wrench_normal, p_center = normal_wrench_patch_center(frameʷ, b)

    IM = MMatrix(I + zeros(SMatrix{4,4,T,16}))
    IM[13:15] .+= SVector{3,T}(p_center.v)
    x_rw_rϕ = Transform3D(FRAME_ϕ, frameʷ, IM)
    x_r2_rϕ = b.x_r²_rʷ * x_rw_rϕ
    x_rϕ_r2 = inv(x_r2_rϕ)

    wrench_stick_phi = transform(wrench_stick_2, x_rϕ_r2)
    wrench_bristle_phi = transform(wrench²_un, x_rϕ_r2)
    wrench_corr_phi = stiction_promoting_soft_clamp(BF.fric_pro, wrench_stick_phi, wrench_bristle_phi)
    wrench² = transform(wrench_corr_phi, x_r2_rϕ)

    ######################

    # NOTE:
    # δδ = -τ (δ + K⁻¹ f²)
    # δδ = -τ (δ + U⁻¹ U⁻ᵀ f²)  # recall: K⁻¹ = U⁻¹ U⁻ᵀ
    # ΔΔ = -τ U (δ + U⁻¹ U⁻ᵀ f²)  # recall: ΔΔ = U δδ
    # ΔΔ = -τ (U δ + U⁻ᵀ f²)
    # ΔΔ = -τ (Δ + U⁻ᵀ f²)
    ΔΔ = get_bristle_d1(tm, bristle_id)
    ΔΔ .= -BF.τ * (Δ  + U⁻¹' * as_static_vector(wrench²))

    return wrench_normal + transform(wrench², b.x_rʷ_r²)
end

spatial_vel_formula(v::SVector{6,T}, b::SVector{3,T}) where {T} = last_3_of_6(v) + cross(first_3_of_6(v), b)

function calc_spatial_bristle_force(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions,
        frameᶜ::CartesianFrame3D, twist_r¹_r²_r²::Twist{T}, δ²::SVector{6,T}) where {N,T}

    # TODO: do this in a frame rotated with the world frame and coincident with the body frame

    b = tm.bodyBodyCache
    tc = b.TractionCache
    BF = c_ins.BristleFriction
    τ⁻¹ = 1 / BF.τ
    μ = c_ins.μ_pair
    frameʷ = tm.frame_world
    k̄ = BF.k̄
    vʳᵉˡ = as_static_vector(twist_r¹_r²_r²)
    wrench_sum = zero(Wrench{T}, frameᶜ)

    for k = 1:length(tc)
        trac = tc.vec[k]
        n̂ = trac.n̂
        n̂² = transform(trac.n̂, b.x_r²_rʷ)
        for k_qp = 1:N
            r² = transform(trac.r_cart[k_qp], b.x_r²_rʷ)
            x̄_δ² = spatial_vel_formula(δ², r².v)
            x̄x̄_vʳᵉˡ = spatial_vel_formula(vʳᵉˡ, r².v)
            term = x̄_δ² + τ⁻¹ * x̄x̄_vʳᵉˡ
            p_dA = calc_p_dA(trac, k_qp)
            λ_s = -k̄ * p_dA * term
            λ_s = vec_sub_vec_proj(λ_s, n̂².v)
            norm_λ_s = norm(λ_s)
            max_fric = μ * p_dA
            if max_fric < norm_λ_s
                λ_s = λ_s * (max_fric / norm_λ_s)
            end
            wrench_sum += Wrench(r², FreeVector3D(frameᶜ, λ_s))
        end
    end
    return wrench_sum
end

function stiction_promoting_soft_clamp(fric_pro::Float64, w_stick::Wrench{T}, w_bristle::Wrench{T}) where {T}

    @framecheck(w_stick.frame, w_bristle.frame)
    corrected_ang = smooth_c1_ramp.(fric_pro * angular(w_bristle), angular(w_stick))
    return Wrench{T}(w_bristle.frame, corrected_ang, linear(w_bristle))
end
