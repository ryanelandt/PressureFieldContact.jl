
function yes_contact!(fric_type::Regularized, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    frame = tm.frame_world
    wrench = zero(Wrench{T}, frame)
    b = tm.bodyBodyCache
    v_tol⁻¹ = c_ins.FrictionModel.v_tol⁻¹
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            cart_vel = trac.v_cart[k]
            n̂ = trac.n̂
            cart_vel_t = vec_sub_vec_proj(cart_vel, n̂)
            mag_vel_t = safe_norm(cart_vel_t.v)
            μ_reg = b.μ * fastSigmoid(mag_vel_t, v_tol⁻¹)
            p_dA = calc_p_dA(trac, k)
            traction_k = p_dA * (n̂ - μ_reg * safe_normalize(cart_vel_t))
            wrench += Wrench(trac.r_cart[k], traction_k)
        end
    end
    return wrench
end

##########################

no_contact!(fric_type::Regularized, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T} = nothing
function no_contact!(fric_type::Bristle, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    BF = c_ins.FrictionModel
    bristle_id = BF.BristleID
    get_bristle_d1(tm, bristle_id) .= -(1 / BF.τ) * get_bristle_d0(tm, bristle_id)
end

function decompose_stiffness!(s::spatialStiffness)
    trace_11 = sqrt(0.33333333333333333 * sum(diag(s.K)[1:3]))
    trace_22 = sqrt(0.33333333333333333 * sum(diag(s.K)[4:6]))
    s.S = Diagonal(SVector(trace_11, trace_11, trace_11, trace_22, trace_22, trace_22))
    s.S⁻¹ = inv(s.S)
    K̄ = s.S⁻¹ * s.K * s.S⁻¹
    mean_val = sum(diag(K̄)) / 6
    s.K̄ = Hermitian(K̄ + I * mean_val * 1.0e-8)
    s.Ū = cholesky(s.K̄).U
    s.Ū⁻¹ = inv(s.Ū)
end

calc_δ(s::spatialStiffness{T}, Δ) where {T} = s.S⁻¹ * (s.Ū⁻¹ * Δ)
calc_Δ(s::spatialStiffness{T}, δ) where {T} = s.Ū * (s.S * δ)

function calc_ΔΔ(τ::Float64, Δ::SVector{6,T}, s::spatialStiffness{T}, wrench_ϕ::Wrench{T}) where {T}
    wrench_ϕ = as_static_vector(wrench_ϕ)
    τ⁻¹ = 1.0 / τ
    return -τ⁻¹ * (Δ + s.Ū⁻¹ * (s.S⁻¹ * wrench_ϕ))
end

function yes_contact!(fric_type::Bristle, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    wrenchʷ_normal, wrenchʷ_fric = bristle_wrench_in_world(tm, c_ins)
    return wrenchʷ_normal + wrenchʷ_fric
end

function transform_δ(v::SVector{6,T}, x) where {T}
    ang = SVector{3,T}(v[1], v[2], v[3])
    lin = SVector{3,T}(v[4], v[5], v[6])
    ang, lin = RigidBodyDynamics.Spatial.transform_spatial_motion(ang, lin, rotation(x), translation(x))
    return vcat(ang, lin)
end

function bristle_wrench_in_world(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    # TODO: improve this to resist slipping more "accurately"
    # TODO: make work for 2 compliant objects. Currently 1 is rigid and 2 is compliant

    BF = c_ins.FrictionModel
    bristle_id = BF.BristleID
    b = tm.bodyBodyCache
    s = b.spatialStiffness

    wrenchʷ_normal, pʷ_center = normal_wrench_patch_center(b)
    calc_patch_spatial_stiffness!(tm, BF)
    transform_stiffness!(s, b.x_r²_rʷ)

    decompose_stiffness!(s)
    Δ² = get_bristle_d0(tm, bristle_id)
    δ² = calc_δ(s, Δ²)
    δʷ = transform_δ(δ², b.x_rʷ_r²)
    wrenchʷ_fric = calc_spatial_bristle_force(tm, c_ins, δʷ, b.twist_r¹_r²)
    wrench²_fric = transform(wrenchʷ_fric, b.x_r²_rʷ)
    get_bristle_d1(tm, bristle_id) .= calc_ΔΔ(BF.τ, Δ², s, wrench²_fric)
    return wrenchʷ_normal, wrenchʷ_fric
end

function set_Δ_from_δʷ(mech_scen, c_ins::ContactInstructions, δʷ)
    SoftContact.force_single_elastic_intersection!(mech_scen, mech_scen.float, c_ins)
    δ² = SoftContact.transform_δ(δʷ, mech_scen.float.bodyBodyCache.x_rʷ_r²)
    Δ = SoftContact.calc_Δ(mech_scen.float.bodyBodyCache.spatialStiffness, δ²)
    b_id = c_ins.FrictionModel.BristleID
    mech_scen.float.s.segments[b_id] .= Δ
end

spatial_vel_formula(v::SVector{6,T}, b::SVector{3,T}) where {T} = last_3_of_6(v) + cross(first_3_of_6(v), b)

# function transform_stiffness!(s::spatialStiffness{T}, x_r²_rʷ) where {T}
#     # see table 2.5 in RBDA
#     R = rotation(x_r²_rʷ)
#     rx = vector_to_skew_symmetric(translation(x_r²_rʷ))
#     X_L = zeros(T,6,6)
#     X_L[1:3, 1:3] = R
#     X_L[4:6, 4:6] = R
#     X_R = X_L'
#     X_L[1:3, 4:6] = -R * rx
#     X_R[4:6, 1:3] = rx * R'
#     s.K = X_L * s.Kʷ * X_R
# end

function transform_stiffness!(s::spatialStiffness{T}, xform) where {T}
    R = rotation(xform)
    t = translation(xform)
    rx = vector_to_skew_symmetric(t)
    X = vcat(hcat(R, (rx * R)), hcat(zeros(SMatrix{3,3,T,9}), R))
    s.K = X * s.Kʷ * X'
end

function calc_patch_spatial_stiffness!(tm::TypedMechanismScenario{N,T}, BF) where {N,T}
    b = tm.bodyBodyCache
    tc = b.TractionCache
    K_11_sum = zeros(SMatrix{3,3,T,9})
    K_12_sum = zeros(SMatrix{3,3,T,9})
    K_22_sum = zeros(SMatrix{3,3,T,9})
    for k = 1:length(tc)
        trac = tc.vec[k]
        n̂ = trac.n̂
        I_minus_n̂n̂ = I - n̂.v * n̂.v'  # suprisingly fast
        for k_qp = 1:N
            p_dA = calc_p_dA(trac, k_qp)
            r = trac.r_cart[k_qp]
            r_skew = vector_to_skew_symmetric(r.v)
            rx_I_minus_n̂n̂ = r_skew * I_minus_n̂n̂
            K_11_sum += -p_dA * rx_I_minus_n̂n̂ * r_skew
            K_12_sum +=  p_dA * rx_I_minus_n̂n̂
            K_22_sum +=  p_dA *    I_minus_n̂n̂
        end
    end
    b.spatialStiffness.Kʷ = BF.k̄ * vcat(hcat(K_11_sum, K_12_sum), hcat(K_12_sum', K_22_sum))
end

function calc_spatial_bristle_force(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions, δ::SVector{6,T},
    twist) where {N,T}

    frame_now = tm.frame_world
    b = tm.bodyBodyCache
    tc = b.TractionCache
    BF = c_ins.FrictionModel
    τ = BF.τ
    μ = b.μ
    k̄ = BF.k̄
    vʳᵉˡ = as_static_vector(twist)
    wrench_sum = zero(Wrench{T}, frame_now)
    for k = 1:length(tc)
        trac = tc.vec[k]
        n̂ = trac.n̂
        for k_qp = 1:N
            r = trac.r_cart[k_qp]
            x̄_δ = spatial_vel_formula(δ, r.v)
            x̄x̄_vʳᵉˡ = spatial_vel_formula(vʳᵉˡ, r.v)
            p_dA = calc_p_dA(trac, k_qp)
            λ_s = -k̄ * p_dA * (x̄_δ + τ * x̄x̄_vʳᵉˡ)
            λ_s = vec_sub_vec_proj(λ_s, n̂.v)
            norm_λ_s = norm(λ_s)
            max_fric = μ * p_dA
            if max_fric < norm_λ_s
                λ_s = λ_s * (max_fric / norm_λ_s)
            end
            wrench_sum += Wrench(r, FreeVector3D(frame_now, λ_s))
        end
    end

    return wrench_sum
end


# function calc_patch_coordinate_system(b::TypedElasticBodyBodyCache{N,T}) where {N,T}
#     frameᶜ = b.mesh_2.FrameID
#     wrenchʷ_normal, pʷ_center = normal_wrench_patch_center(b)
#     p_centerᶜ = transform(pʷ_center, b.x_r²_rʷ)
#     ### want to calculate in a frame aligned with the compliant body frame and instaneously coincident with the COP
#     IM = I + zeros(MMatrix{4,4,T,16})
#     IM[13:15] .+= SVector{3,T}(-p_centerᶜ.v)
#     x_rϕ_rʷ = Transform3D(frameᶜ, FRAME_ϕ, IM) * b.x_r²_rʷ
#     return x_rϕ_rʷ, wrenchʷ_normal
# end

# function stiction_promoting_soft_clamp(fric_pro::Float64, w_stick::Wrench{T}, w_bristle::Wrench{T}) where {T}
#     @framecheck(w_stick.frame, w_bristle.frame)
#     corrected_ang = smooth_c1_ramp.(fric_pro * angular(w_bristle), angular(w_stick))
#     return Wrench{T}(w_bristle.frame, corrected_ang, linear(w_bristle))
# end

# K̇_11_sum = zeros(SMatrix{3,3,T,9})
# K̇_12_sum = zeros(SMatrix{3,3,T,9})
# K̇_22_sum = zeros(SMatrix{3,3,T,9})
#
# DOT_n̂² = cross(angular(twist_cf), n̂².v)  # works because r2 is rigid  works because point is treated as fixed at this instant
# DOT_r²_skew = vector_to_skew_symmetric(point_velocity(twist_cf, r²).v)
# DOT_I_minus_n̂n̂ = -(DOT_n̂² * n̂².v' + n̂².v * DOT_n̂²')
# DOT_rx_I_minus_n̂n̂ = DOT_I_minus_n̂n̂ * r²_skew + I_minus_n̂n̂ * DOT_r²_skew
# K̇_11_sum += -p_dA * (DOT_rx_I_minus_n̂n̂ * r²_skew + rx_I_minus_n̂n̂ * DOT_r²_skew)
# K̇_12_sum +=  p_dA * DOT_rx_I_minus_n̂n̂
# K̇_22_sum +=  p_dA * DOT_I_minus_n̂n̂
#
# K̇ = (BF.k̄ * μ) * vcat(hcat(K̇_11_sum, K̇_12_sum), hcat(K̇_12_sum', K̇_22_sum))

# δδϕ_star = -BF.τ * (δϕ + sqrt_K⁻¹ * (sqrt_K⁻¹ * as_static_vector(wrench_ϕ)) )
# γγ .= sqrt_K * δδϕ_star

# b.K.data .= K
# K⁻¹ = inv(K)
# cf = cholesky(b.K)
# U = cf.U
#
# # U = sqrtm(K)
# U⁻¹ = inv(U)
#
# Δϕ = get_bristle_d0(tm, bristle_id)
# Δϕ = SVector{6,T}(Δϕ[1], Δϕ[2], Δϕ[3], Δϕ[4], Δϕ[5], Δϕ[6])
# δϕ = U⁻¹ * Δϕ
# δϕ = SVector{6,T}(δϕ[1], δϕ[2], δϕ[3], δϕ[4], δϕ[5], δϕ[6])
#
#
#
# ΔΔϕ = get_bristle_d1(tm, bristle_id)
# ΔΔϕ .= U * δδϕ_star

# b.K.data .= K
# K⁻¹ = inv(K)
# cf = cholesky(b.K)
# U = cf.U
#
# # U = sqrtm(K)
# U⁻¹ = inv(U)
#
# Δϕ = get_bristle_d0(tm, bristle_id)
# Δϕ = SVector{6,T}(Δϕ[1], Δϕ[2], Δϕ[3], Δϕ[4], Δϕ[5], Δϕ[6])
# δϕ = U⁻¹ * Δϕ
# δϕ = SVector{6,T}(δϕ[1], δϕ[2], δϕ[3], δϕ[4], δϕ[5], δϕ[6])
#
#
# wrench_ϕ = calc_spatial_bristle_force_cf(tm, c_ins, δϕ, twist_r¹_r²_rϕ, x_rϕ_rʷ)
# δδϕ_star = -BF.τ * (δϕ + K⁻¹ * as_static_vector(wrench_ϕ))
#
# ΔΔϕ = get_bristle_d1(tm, bristle_id)
# ΔΔϕ .= U * δδϕ_star

#     U⁻¹ = SMatrix{6,6,T,36}(inv(U))
#     δ = U⁻¹ * Δ

# K⁻¹ = inv(K)
#
#
# τϕ = τϕ * bbb
# # τ = K δ
# # δ = K⁻¹ τ
# δϕ =  K⁻¹ * τϕ
#
#
# δδϕ_work = -0.5 * K⁻¹ * K̇ * δϕ
# δδϕ = δδϕ_star + δδϕ_work
# # ττϕ .= K * (δδϕ_star + δδϕ_work)
# ττϕ = get_bristle_d1(tm, bristle_id)
# ττϕ .= K * δδϕ + K̇ * δϕ
# ττϕ .= ττϕ * (1 / bbb)

# bbb = 1.0e6
#
# τϕ = get_bristle_d0(tm, bristle_id)
# τϕ = SVector{6,T}(τϕ[1], τϕ[2], τϕ[3], τϕ[4], τϕ[5], τϕ[6])
# τϕ = τϕ * bbb
# # τ = K δ
# # δ = K⁻¹ τ
# δϕ =  K⁻¹ * τϕ
#
# wrench_ϕ = calc_spatial_bristle_force_cf(tm, c_ins, δϕ, twist_r¹_r²_rϕ, x_rϕ_rʷ)
#
# δδϕ_work = -0.5 * K⁻¹ * K̇ * δϕ
# δδϕ_star = -BF.τ * (δϕ + K⁻¹ * as_static_vector(wrench_ϕ))
# δδϕ = δδϕ_star + δδϕ_work
# # ττϕ .= K * (δδϕ_star + δδϕ_work)
# ττϕ = get_bristle_d1(tm, bristle_id)
# ττϕ .= K * δδϕ + K̇ * δϕ
# ττϕ .= ττϕ * (1 / bbb)

# function calc_spatial_bristle_force_world(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions, δʷ::SVector{6,T},
#     twist_r¹_r²_rʷ::Twist{T}, ) where {N,T}
#
#     b = tm.bodyBodyCache
#     tc = b.TractionCache
#     BF = c_ins.FrictionModel
#     τ⁻¹ = 1 / BF.τ
#     μ = b.μ
#     frameʷ = tm.frame_world
#     k̄ = BF.k̄
#     vʳᵉˡ = as_static_vector(twist_r¹_r²_rʷ)
#     wrench_sum = zero(Wrench{T}, frameʷ)
#     for k = 1:length(tc)
#         trac = tc.vec[k]
#         n̂ʷ = trac.n̂
#         for k_qp = 1:N
#             rʷ = trac.r_cart[k_qp]
#             x̄_δʷ = spatial_vel_formula(δʷ, rʷ.v)
#             x̄x̄_vʳᵉˡ = spatial_vel_formula(vʳᵉˡ, rʷ.v)
#             term = x̄_δʷ + τ⁻¹ * x̄x̄_vʳᵉˡ
#             p_dA = calc_p_dA(trac, k_qp)
#             λ_s = -k̄ * p_dA * term
#             λ_s = vec_sub_vec_proj(λ_s, n̂ʷ.v)
#             norm_λ_s = norm(λ_s)
#             max_fric = μ * p_dA
#             if max_fric < norm_λ_s
#                 λ_s = λ_s * (max_fric / norm_λ_s)
#             end
#             wrench_sum += Wrench(rʷ, FreeVector3D(frameʷ, λ_s))
#         end
#     end
#     wrench_sum = transform(wrench_sum, b.x_r²_rʷ)
#     return wrench_sum
# end

# function calc_spatial_bristle_force(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions,
#     frameᶜ::CartesianFrame3D, twist_r¹_r²_r²::Twist{T}, δ²::SVector{6,T}) where {N,T}
#
#     # TODO: do this in a frame rotated with the world frame and coincident with the body frame
#
#     b = tm.bodyBodyCache
#     tc = b.TractionCache
#     BF = c_ins.FrictionModel
#     τ⁻¹ = 1 / BF.τ
#     μ = b.μ
#     frameʷ = tm.frame_world
#     k̄ = BF.k̄
#     vʳᵉˡ = as_static_vector(twist_r¹_r²_r²)
#     wrench_sum = zero(Wrench{T}, frameᶜ)
#
#     for k = 1:length(tc)
#         trac = tc.vec[k]
#         n̂ = trac.n̂
#         n̂² = transform(trac.n̂, b.x_r²_rʷ)
#         for k_qp = 1:N
#             r² = transform(trac.r_cart[k_qp], b.x_r²_rʷ)
#             x̄_δ² = spatial_vel_formula(δ², r².v)
#             x̄x̄_vʳᵉˡ = spatial_vel_formula(vʳᵉˡ, r².v)
#             term = x̄_δ² + τ⁻¹ * x̄x̄_vʳᵉˡ
#             p_dA = calc_p_dA(trac, k_qp)
#             λ_s = -k̄ * p_dA * term
#             λ_s = vec_sub_vec_proj(λ_s, n̂².v)
#             norm_λ_s = norm(λ_s)
#             max_fric = μ * p_dA
#             if max_fric < norm_λ_s
#                 λ_s = λ_s * (max_fric / norm_λ_s)
#             end
#             wrench_sum += Wrench(r², FreeVector3D(frameᶜ, λ_s))
#         end
#     end
#     return wrench_sum
# end

# function calc_patch_spatial_stiffness_and_derivative(tm::TypedMechanismScenario{N,T}, BF, c_ins::ContactInstructions, twist_r¹_r²_r²) where {N,T}
#     # TODO: make work for 2 compliant objects. Currently 1 is rigid and 2 is compliant
#     # TODO: do this in a frame rotated with the world frame and coincident with the body frame
#     # [-rx I_minus_n̂n̂ rx,  rx I_minus_n̂n̂;
#     #                         I_minus_n̂n̂]
#
#     b = tm.bodyBodyCache
#     frameᶜ = b.mesh_2.FrameID
#     tc = b.TractionCache
#     K_11_sum = zeros(SMatrix{3,3,T,9})
#     K_12_sum = zeros(SMatrix{3,3,T,9})
#     K_22_sum = zeros(SMatrix{3,3,T,9})
#     K̇_11_sum = zeros(SMatrix{3,3,T,9})
#     K̇_12_sum = zeros(SMatrix{3,3,T,9})
#     K̇_22_sum = zeros(SMatrix{3,3,T,9})
#     for k = 1:length(tc)
#         trac = tc.vec[k]
#         n̂² = transform(trac.n̂, b.x_r²_rʷ)
#         I_minus_n̂n̂ = I - n̂².v * n̂².v'  # suprisingly fast
#         for k_qp = 1:N
#             p_dA = calc_p_dA(trac, k_qp)
#             r² = transform(trac.r_cart[k_qp], b.x_r²_rʷ)
#             r²_skew = vector_to_skew_symmetric(r².v)
#             rx_I_minus_n̂n̂ = r²_skew * I_minus_n̂n̂
#             K_11_sum += -p_dA * rx_I_minus_n̂n̂ * r²_skew
#             K_12_sum +=  p_dA * rx_I_minus_n̂n̂
#             K_22_sum +=  p_dA *    I_minus_n̂n̂
#             DOT_n̂² = cross(angular(twist_r¹_r²_r²), n̂².v)  # works because r2 is rigid
#             # works because point is treated as fixed at this instant
#             DOT_r²_skew = vector_to_skew_symmetric(point_velocity(twist_r¹_r²_r², r²).v)
#             DOT_I_minus_n̂n̂ = -(DOT_n̂² * n̂².v' + n̂².v * DOT_n̂²')
#             DOT_rx_I_minus_n̂n̂ = DOT_I_minus_n̂n̂ * r²_skew + I_minus_n̂n̂ * DOT_r²_skew
#             K̇_11_sum += -p_dA * (DOT_rx_I_minus_n̂n̂ * r²_skew + rx_I_minus_n̂n̂ * DOT_r²_skew)
#             K̇_12_sum +=  p_dA * DOT_rx_I_minus_n̂n̂
#             K̇_22_sum +=  p_dA * DOT_I_minus_n̂n̂
#         end
#     end
#     μ = b.μ
#     K = (BF.k̄ * μ) * vcat(hcat(K_11_sum, K_12_sum), hcat(K_12_sum', K_22_sum))
#     K̇ = (BF.k̄ * μ) * vcat(hcat(K̇_11_sum, K̇_12_sum), hcat(K̇_12_sum', K̇_22_sum))
#     # TODO: choose this in a principled way
#     K += SMatrix{6,6,T,36}(diagm(0=>[1.0e-6, 1.0e-6, 1.0e-6, 1.0e-3, 1.0e-3, 1.0e-3]))
#     return K, K̇
# end

# function calc_patch_spatial_stiffness(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions, K) where {N,T}
#
#     # TODO: do this in a frame rotated with the world frame and coincident with the body frame
#
#     b = tm.bodyBodyCache
#     frameᶜ = b.mesh_2.FrameID
#     tc = b.TractionCache
#
#     Q_11_sum = zeros(SMatrix{3,3,T,9})
#     Q_12_sum = zeros(SMatrix{3,3,T,9})
#     Q_22_sum = zeros(SMatrix{3,3,T,9})
#     for k = 1:length(tc)
#         trac = tc.vec[k]
#         n̂² = transform(trac.n̂, b.x_r²_rʷ)
#         Q_22 = I - n̂².v * n̂².v'  # suprisingly fast
#         for k_qp = 1:N
#             r² = transform(trac.r_cart[k_qp], b.x_r²_rʷ)
#             r²_skew = vector_to_skew_symmetric(r².v)
#             Q_12 = r²_skew * Q_22
#             minus_Q_11 = Q_12 * r²_skew
#             p_dA = calc_p_dA(trac, k_qp)
#             Q_11_sum -= p_dA * minus_Q_11
#             Q_12_sum += p_dA * Q_12
#             Q_22_sum += p_dA * Q_22
#         end
#     end
#     μ = b.μ
#     K[1:3, 1:3] .= μ * Q_11_sum
#     K[1:3, 4:6] .= μ * Q_12_sum
#     K[4:6, 4:6] .= μ * Q_22_sum
#     return nothing
# end

# function veil_friction!(frameʷ::CartesianFrame3D, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
#     # TODO: do all calculations relative to patch center
#
#     BF = c_ins.FrictionModel
#     bristle_id = BF.BristleID
#     Δ = get_bristle_d0(tm, bristle_id)
#     Δ = SVector{6,T}(Δ[1], Δ[2], Δ[3], Δ[4], Δ[5], Δ[6])
#
#     b = tm.bodyBodyCache
#     twist_r¹_r² = b.twist_r¹_r²
#     twist_r¹_r²_r² = transform(twist_r¹_r², b.x_r²_rʷ)
#
#     frameᶜ = b.mesh_2.FrameID
#     @framecheck(twist_r¹_r²_r².frame, frameᶜ)
#     v_spatial_rel = as_static_vector(twist_r¹_r²_r²)
#
#     K = b.K
#     fill!(K.data, zero(T))
#     calc_patch_spatial_stiffness(tm, c_ins, K.data)
#
#     make_stiffness_PD!(K.data)
#     K.data .*= BF.k̄
#
#     K_err = Matrix(K)
#     if isposdef(K_err)
#         #
#     else
#         println("")
#         for k_33 = 1:6
#             println(K_err[k_33, 1:6])
#         end
#         println("")
#         println("smoking_gun = [")
#         for k = 1:36
#             a4 = K[k]
#             @printf("%.17f", a4)
#             (k == 36) || (print(","))
#         end
#         println("]")
#         println("")
#         println("length(b.TractionCache): ", length(b.TractionCache))
#         println("")
#         error("something is very very wrong")
#     end
#
#     # NOTE:
#     # Forward       Inverse
#     # K = Uᵀ U      K⁻¹ = U⁻¹ U⁻ᵀ
#     # Δ = U δ       δ = U⁻¹ Δ
#     # ΔΔ = U δδ     δδ = U⁻¹ ΔΔ
#
#     cf = cholesky(K)
#     U = cf.U
#     U⁻¹ = SMatrix{6,6,T,36}(inv(U))
#     δ = U⁻¹ * Δ
#     wrench²_un = calc_spatial_bristle_force(tm, c_ins, frameᶜ, twist_r¹_r²_r², δ)
#
#     ######################
#
#     τ_s2 = -K * (Δ + BF.τ * v_spatial_rel)
#     wrench_stick_2 = Wrench{T}(frameᶜ, SVector{3,T}(τ_s2[1], τ_s2[2], τ_s2[3]), SVector{3,T}(τ_s2[4], τ_s2[5], τ_s2[6]))
#     wrench_normal, p_center = normal_wrench_patch_center(frameʷ, b)
#
#     IM = MMatrix(I + zeros(SMatrix{4,4,T,16}))
#     IM[13:15] .+= SVector{3,T}(p_center.v)
#     x_rw_rϕ = Transform3D(FRAME_ϕ, frameʷ, IM)
#     x_r2_rϕ = b.x_r²_rʷ * x_rw_rϕ
#     x_rϕ_r2 = inv(x_r2_rϕ)
#
#     wrench_stick_phi = transform(wrench_stick_2, x_rϕ_r2)
#     wrench_bristle_phi = transform(wrench²_un, x_rϕ_r2)
#     wrench_corr_phi = stiction_promoting_soft_clamp(BF.fric_pro, wrench_stick_phi, wrench_bristle_phi)
#     wrench² = transform(wrench_corr_phi, x_r2_rϕ)
#
#     ######################
#
#     # NOTE:
#     # δδ = -τ (δ + K⁻¹ f²)
#     # δδ = -τ (δ + U⁻¹ U⁻ᵀ f²)  # recall: K⁻¹ = U⁻¹ U⁻ᵀ
#     # ΔΔ = -τ U (δ + U⁻¹ U⁻ᵀ f²)  # recall: ΔΔ = U δδ
#     # ΔΔ = -τ (U δ + U⁻ᵀ f²)
#     # ΔΔ = -τ (Δ + U⁻ᵀ f²)
#     ΔΔ = get_bristle_d1(tm, bristle_id)
#     ΔΔ .= -BF.τ * (Δ  + U⁻¹' * as_static_vector(wrench²))
#
#     return wrench_normal + transform(wrench², b.x_rʷ_r²)
# end
