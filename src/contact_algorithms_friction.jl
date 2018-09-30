function regularized_friction(frame::CartesianFrame3D, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    wrench = zeroWrench(frame, T)
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            cart_vel_crw_t = trac.v_cart_t[k]
            mag_vel_t  = safeNorm(cart_vel_crw_t.v)
            mu_reg     = b.mu * fastSigmoid(mag_vel_t)
            traction_k = trac.p_dA[k] * (trac.traction_normal - mu_reg * safeNormalize(cart_vel_crw_t))
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
