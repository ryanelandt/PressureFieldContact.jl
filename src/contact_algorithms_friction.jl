function friction_regularization end
function friction_bristle end

function friction_model(::typeof(friction_regularization), wrench::Wrench{T}, b::TypedElasticBodyBodyCache{N,T}) where {N,T}
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
