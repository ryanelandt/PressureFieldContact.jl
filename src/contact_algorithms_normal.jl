
function normal_wrench(b::TypedElasticBodyBodyCache{N,T}) where {N,T}
    wrench_lin = zeros(SVector{3,T})
    wrench_ang = zeros(SVector{3,T})
    @inbounds begin
    for k_trac = 1:length(b.TractionCache)
        trac = b.TractionCache[k_trac]
        # for k = 1:N
            p_dA = calc_p_dA(trac)
            λ_s = -p_dA * trac.n̂.v
            wrench_lin += λ_s
            wrench_ang += cross(trac.r_cart.v, λ_s)
        # end
    end
    end
    return Wrench(b.mesh_2.FrameID, wrench_ang, wrench_lin)
end
