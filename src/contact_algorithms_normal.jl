
function normal_wrench(b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    frame = b.frame_world
    wrench = zero(Wrench{T}, frame)
    @inbounds begin
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            p_dA = calc_p_dA(trac, k)
            wrench += Wrench(trac.r_cart[k], p_dA * trac.n̂)
        end
    end
    end
    return wrench
end

function normal_wrench_patch_center(b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    frame = b.frame_world
    wrench = zero(Wrench{T}, frame)
    int_p_r_dA = zeros(SVector{3,T})
    int_p_dA = zero(T)
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        for k = 1:N
            p_dA = calc_p_dA(trac, k)
            wrench += Wrench(trac.r_cart[k], p_dA * trac.n̂)
            int_p_r_dA += trac.r_cart[k].v * p_dA
            int_p_dA += p_dA
        end
    end
    p_center = Point3D(frame, int_p_r_dA / int_p_dA)
    return wrench, p_center
end
