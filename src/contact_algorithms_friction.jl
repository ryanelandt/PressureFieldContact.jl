
function yes_contact!(fric_type::Regularized, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    b = tm.bodyBodyCache
    frame = b.mesh_2.FrameID
    v_tol⁻¹ = c_ins.FrictionModel.v_tol⁻¹
    v_tol = 1 / v_tol⁻¹
    wrench_lin = zeros(SVector{3,T})
    wrench_ang = zeros(SVector{3,T})
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        cart_vel = trac.v_cart
        n̂ = trac.n̂
        cart_vel_t = vec_sub_vec_proj(cart_vel, n̂)
        p_dA = calc_p_dA(trac)

        mag_vel_t² = dot(cart_vel_t, cart_vel_t)
        if mag_vel_t² <= v_tol^2
            term = b.μ * cart_vel_t * v_tol⁻¹
        else
            term = b.μ * cart_vel_t / sqrt(mag_vel_t²)
        end

        traction_k = p_dA * (term - n̂)
        wrench_lin += traction_k
        wrench_ang += cross(trac.r_cart, traction_k)
    end
    return Wrench(frame, wrench_ang, wrench_lin)
end

##########################

no_contact!(fric_type::Regularized, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T} = nothing
function no_contact!(fric_type::Bristle, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
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

function bristle_wrench_in_world(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    BF = c_ins.FrictionModel
    bristle_id = BF.BristleID
    b = tm.bodyBodyCache
    spatial_stiffness = b.spatialStiffness
    τ⁻¹ = 1 / BF.τ
    calc_patch_spatial_stiffness!(tm, BF)
    calc_K_sqrt⁻¹!(spatial_stiffness)
    s = get_bristle_d0(tm, bristle_id)
    Δ² = spatial_stiffness.K⁻¹_sqrt * s
    wrench²_normal, wrench²_fric = calc_spatial_bristle_force(tm, c_ins, Δ², b.twist_r²_r¹_r²)
    f_spatial = as_static_vector(wrench²_fric)
    get_bristle_d1(tm, bristle_id) .= τ⁻¹ * ( spatial_stiffness.K⁻¹_sqrt * f_spatial - s)
    return wrench²_normal, -wrench²_fric
end

#########################################################

function yes_contact!(fric_type::Bristle, tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions) where {N,T}
    wrench²_normal, wrench²_fric = bristle_wrench_in_world(tm, c_ins)
    return wrench²_normal + wrench²_fric
end

spatial_vel_formula(v::SVector{6,T}, b::SVector{3,T}) where {T} = last_3_of_6(v) + cross(first_3_of_6(v), b)

function calc_patch_spatial_stiffness!(tm::TypedMechanismScenario{N,T}, BF) where {N,T}
    b = tm.bodyBodyCache
    tc = b.TractionCache
    K_11_sum = zeros(SMatrix{3,3,T,9})
    K_12_sum = zeros(SMatrix{3,3,T,9})
    K_22_sum = zeros(SMatrix{3,3,T,9})
    for k = 1:length(tc)
        trac = tc.vec[k]
        n̂ = trac.n̂
        p_dA = calc_p_dA(trac)
        r = trac.r_cart
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

function calc_spatial_bristle_force(tm::TypedMechanismScenario{N,T}, c_ins::ContactInstructions, δ::SVector{6,T},
    twist) where {N,T}

    # TODO: have this include normal tractions too
    b = tm.bodyBodyCache
    frame_now = b.mesh_2.FrameID
    tc = b.TractionCache
    BF = c_ins.FrictionModel
    τ = BF.τ
    μ = b.μ
    k̄ = BF.k̄
    vʳᵉˡ = as_static_vector(twist)
    wrench_lin = zeros(SVector{3,T})
    wrench_ang = zeros(SVector{3,T})
    normal_wrench_lin = zeros(SVector{3,T})
    normal_wrench_ang = zeros(SVector{3,T})
    for k = 1:length(tc)
        trac = tc.vec[k]
        n̂ = trac.n̂
        r = trac.r_cart
        x̄_δ = spatial_vel_formula(δ, r)
        x̄x̄_vʳᵉˡ = spatial_vel_formula(vʳᵉˡ, r)
        p_dA = calc_p_dA(trac)
        λ_s = k̄ * p_dA * (x̄_δ - τ * x̄x̄_vʳᵉˡ)
        λ_s = vec_sub_vec_proj(λ_s, n̂)
        mag_λ_s = norm(λ_s)
        if μ * p_dA < mag_λ_s
            λ_s = μ * p_dA * λ_s / mag_λ_s
        end
        wrench_lin += λ_s
        wrench_ang += cross(r, λ_s)
        # Normal Wrench
        λₙ = -p_dA * n̂
        normal_wrench_lin += λₙ
        normal_wrench_ang += cross(r, λₙ)
    end
    wrench_fric   = Wrench(frame_now, wrench_ang, wrench_lin)
    wrench_normal = Wrench(frame_now, normal_wrench_ang, normal_wrench_lin)
    return wrench_normal, wrench_fric
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
