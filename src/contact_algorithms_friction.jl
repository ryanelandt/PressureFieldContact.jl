
function calc_clamped_piecewise(x::T, x_1::Float64, x_2::Float64, y_1::Float64, y_2::Float64) where {T <: Real}
	# The function varies linearly between y_1 at x_1 to y_2 at x_2.
	# It is assumed that y_2 < y_1.
	# The output is bounded between y_1 and y_2.

	k = (y_2 - y_1) / (x_2 - x_1)
	y = y_1 + (x - x_1) * k
	return clamp(y, y_2, y_1)
end

function calc_μ(ins::Bristle, δ::T, mag_T̄s::T) where {T}
    δ_max = ins.δ_max
	μs    = ins.μs
	μd    = ins.μd
    γ     = calc_clamped_piecewise(δ, δ_max, ins.δ_μ_is_0, 1.0, 0.0)
    return γ * calc_clamped_piecewise(mag_T̄s, μs, ins.T̄s_μ_is_μd, μs, μd)
end

function calc_μ(ins::Regularized, mag_vel_t::T) where {T}
	μs = ins.μs
	μd = ins.μd
	return calc_clamped_piecewise(mag_vel_t, ins.v_tol, ins.v_μ_is_μd, μs, μd)
end

### Traction
function traction(ins::Bristle, δ::T, T̄s::SVector{3,T}, p_dA::T) where {T}
	mag_T̄s = norm(T̄s)
	μ = calc_μ(ins, δ, mag_T̄s)
	if μ == 0
		return 0 * T̄s
	elseif mag_T̄s <= μ
		T̄c = T̄s
	else
		T̄c = T̄s * μ * (1 / mag_T̄s)
	end
	return T̄c * p_dA
end

function traction(ins::Regularized, cart_vel_t::SVector{3,T}, p_dA::T) where {T}
	mag_vel_t = norm(cart_vel_t)
	μ = calc_μ(ins, mag_vel_t)
	T_c = -μ * p_dA * cart_vel_t
	if mag_vel_t <= ins.v_tol
		return T_c * (1 / ins.v_tol)
	else
		return T_c * (1 / mag_vel_t)
	end
end

function yes_contact!(fric_type::Regularized, tm::TypedMechanismScenario{T}, c_ins::ContactInstructions) where {T}
    b = tm.bodyBodyCache
    frame = b.mesh_2.FrameID
    # v_tol⁻¹ = c_ins.FrictionModel.v_tol⁻¹
    # v_tol = 1 / v_tol⁻¹
    twist = b.twist_r²_r¹_r²
    vʳᵉˡ = as_static_vector(twist)
    wrench_lin = zeros(SVector{3,T})
    wrench_ang = zeros(SVector{3,T})
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        r = trac.r_cart
        n̂ = trac.n̂
        cart_vel = spatial_vel_formula(vʳᵉˡ, r)
        cart_vel_t = vec_sub_vec_proj(cart_vel, n̂)
        p_dA = calc_p_dA(trac)

        # T_c = regularized_μs_μd(cart_vel_t, p_dA, v_tol, b.μs, b.μd)

		T_c = traction(fric_type, cart_vel_t, p_dA)

        traction_k = -p_dA * n̂ - T_c

        wrench_lin += traction_k
        wrench_ang += cross(r, traction_k)
    end
    return Wrench(frame, wrench_ang, wrench_lin)
end

##########################

no_contact!(fric_type::Regularized, tm::TypedMechanismScenario{T}, c_ins::ContactInstructions) where {T} = nothing
function no_contact!(fric_type::Bristle, tm::TypedMechanismScenario{T}, c_ins::ContactInstructions) where {T}
    BF = c_ins.FrictionModel
    bristle_id = BF.BristleID
    get_bristle_d1(tm, bristle_id) .= -(1 / BF.τ) * get_bristle_d0(tm, bristle_id)
end

#########################################################

# TODO: make this allocate less
function calc_K_sqrt⁻¹!(s::spatialStiffness{T}) where {T}
    EF = eigen!(s.K)
    s.σ_sqrt.diag .= EF.values
    max_σ = maximum(s.σ_sqrt.diag)
    for k = 1:6
        s.σ_sqrt.diag[k] = 1.0 / sqrt(max(s.σ_sqrt.diag[k], max_σ * 1.0e-16))
    end
    mul!(s.mul_pre, EF.vectors, s.σ_sqrt)
    mul!(s.K⁻¹_sqrt, s.mul_pre, EF.vectors')
end

function bristle_wrench_in_world(tm::TypedMechanismScenario{T}, c_ins::ContactInstructions) where {T}
    # cop: patch center of pressure

    BF = c_ins.FrictionModel
    bristle_id = BF.BristleID
    b = tm.bodyBodyCache
    spatial_stiffness = b.spatialStiffness
    τ⁻¹ = 1 / BF.τ
    cop, wrench²_normal = normal_wrench_cop(b)
    calc_patch_spatial_stiffness!(tm, BF, cop.v)
    calc_K_sqrt⁻¹!(spatial_stiffness)
    s = get_bristle_d0(tm, bristle_id)
    Δ² = spatial_stiffness.K⁻¹_sqrt * s
    wrenchᶜᵒᵖ_fric, wrench²_fric = calc_spatial_bristle_force(tm, c_ins, Δ², b.twist_r²_r¹_r², cop.v)
    get_bristle_d1(tm, bristle_id) .= τ⁻¹ * ( spatial_stiffness.K⁻¹_sqrt * wrenchᶜᵒᵖ_fric - s)
    return wrench²_normal, -wrench²_fric
end

#########################################################

function yes_contact!(fric_type::Bristle, tm::TypedMechanismScenario{T}, c_ins::ContactInstructions) where {T}
    wrench²_normal, wrench²_fric = bristle_wrench_in_world(tm, c_ins)
    return wrench²_normal + wrench²_fric
end

spatial_vel_formula(v::SVector{6,T}, b::SVector{3,T}) where {T} = last_3_of_6(v) + cross(first_3_of_6(v), b)

function calc_patch_spatial_stiffness!(tm::TypedMechanismScenario{T}, BF, cop::SVector{3,T}) where {T}
    b = tm.bodyBodyCache
    tc = b.TractionCache
    K_11_sum = zeros(SMatrix{3,3,T,9})
    K_12_sum = zeros(SMatrix{3,3,T,9})
    K_22_sum = zeros(SMatrix{3,3,T,9})
    for k = 1:length(tc)
        trac = tc.vec[k]
        n̂ = trac.n̂
        p_dA = calc_p_dA(trac)
        r = trac.r_cart - cop
        K_22_sum += p_dA * (I - n̂ * n̂')
        r_x_n̂     = cross(r, n̂)
        K_12_sum += p_dA * (vector_to_skew_symmetric(r) - r_x_n̂ * n̂')
        r_skew²   = vector_to_skew_symmetric_squared(r)
        K_11_sum -= p_dA * (r_skew² + r_x_n̂ * r_x_n̂')
    end
    b.spatialStiffness.K.data[1:3, 1:3] .= K_11_sum
    b.spatialStiffness.K.data[4:6, 1:3] .= K_12_sum'
    b.spatialStiffness.K.data[1:3, 4:6] .= K_12_sum
    b.spatialStiffness.K.data[4:6, 4:6] .= K_22_sum
    b.spatialStiffness.K.data .*= BF.k̄
end

function calc_spatial_bristle_force(tm::TypedMechanismScenario{T}, c_ins::ContactInstructions, Δ²::SVector{6,T},
    twist, cop::SVector{3,T}) where {T}

    b = tm.bodyBodyCache
    tc = b.TractionCache
    BF = c_ins.FrictionModel
    τ = BF.τ
    k̄ = BF.k̄
    vʳᵉˡ = as_static_vector(twist)
    wrench_lin = zeros(SVector{3,T})
    wrench_ang = zeros(SVector{3,T})
    for k = 1:length(tc)
        trac = tc.vec[k]
        n̂ = trac.n̂
        r = trac.r_cart
        x² = r - cop
        δ² = spatial_vel_formula(Δ², x²)
        ṙ²_perp = spatial_vel_formula(vʳᵉˡ, r)
        p_dA = calc_p_dA(trac)

		T̄s = k̄ * (δ² - τ * ṙ²_perp)
        T̄s = vec_sub_vec_proj(T̄s, n̂)
		T_c = traction(BF, norm(δ²), T̄s, p_dA)

        wrench_lin += T_c
        wrench_ang += cross(x², T_c)
    end
    wrenchᶜᵒᵖ_fric = vcat(wrench_ang, wrench_lin)
    wrench²_fric = vcat(wrench_ang + cross(cop, wrench_lin), wrench_lin)
    return wrenchᶜᵒᵖ_fric, Wrench(b.mesh_2.FrameID, wrench_ang + cross(cop, wrench_lin), wrench_lin)
end




# function transform_stiffness!(s::spatialStiffness{T}, xform) where {T}
#     R = rotation(xform)
#     t = translation(xform)
#     rx = vector_to_skew_symmetric(t)
#     X = vcat(hcat(R, (rx * R)), hcat(zeros(SMatrix{3,3,T,9}), R))
#     mul!(s.mul_pre, X, s.K)
#     mul!(s.K.data, s.mul_pre, X')
# end

# function calc_s(s::spatialStiffness, Δ)
#     # this function should not be called in tight loops
#     return inv(s.K⁻¹_sqrt) * Δ
# end

# function transform_δ(v::SVector{6,T}, x) where {T}
#     ang = SVector{3,T}(v[1], v[2], v[3])
#     lin = SVector{3,T}(v[4], v[5], v[6])
#     ang, lin = RigidBodyDynamics.Spatial.transform_spatial_motion(ang, lin, rotation(x), translation(x))
#     return vcat(ang, lin)
# end

# function set_s_from_Δʷ(mech_scen, c_ins::ContactInstructions, δʷ)
#     force_single_elastic_intersection!(mech_scen, mech_scen.float, c_ins)
#     δ² = transform_δ(δʷ, mech_scen.float.bodyBodyCache.x_rʷ_r²)
#     Δ = calc_s(mech_scen.float.bodyBodyCache.spatialStiffness, δ²)
#     b_id = c_ins.FrictionModel.BristleID
#     mech_scen.float.s.segments[b_id] .= Δ
# end
