
function calc_stiffness_finite_difference(r::SVector{3,T}, n̂::SVector{3,T}, tol::Float64=1.0e-8) where {T}
    function make_homogenous(R, t)
        R_t = hcat(R, t)
        bottom = SMatrix{1,4,Float64,4}(0.0, 0.0, 0.0, 1.0)
        return vcat(R_t, bottom)
    end

    function calc_stiffness_element(H::SMatrix{4,4,Float64,16}, r::SVector{3,Float64}, n̂::SVector{3,Float64})
        a1 = H * SVector{4,Float64}(r[1], r[2], r[3], 1.0)
        r_pert = SVector{3,Float64}(a1[1], a1[2], a1[3])
        lin_proj = vec_sub_vec_proj(r_pert - r, n̂)
        ang_proj = cross(r, lin_proj)
        return vcat(ang_proj, lin_proj)
    end

    tol = 1.0e-08
    R_0 = one(RotMatrix{3,Float64})
    sv3_0 = zeros(SVector{3,Float64})
    H_1 = make_homogenous(AngleAxis(tol, 1.0, 0.0, 0.0), sv3_0)
    H_2 = make_homogenous(AngleAxis(tol, 0.0, 1.0, 0.0), sv3_0)
    H_3 = make_homogenous(AngleAxis(tol, 0.0, 0.0, 1.0), sv3_0)
    H_4 = make_homogenous(R_0, SVector{3,Float64}(tol, 0.0, 0.0))
    H_5 = make_homogenous(R_0, SVector{3,Float64}(0.0, tol, 0.0))
    H_6 = make_homogenous(R_0, SVector{3,Float64}(0.0, 0.0, tol))
    K_1 = calc_stiffness_element(H_1, r, n̂)
    K_2 = calc_stiffness_element(H_2, r, n̂)
    K_3 = calc_stiffness_element(H_3, r, n̂)
    K_4 = calc_stiffness_element(H_4, r, n̂)
    K_5 = calc_stiffness_element(H_5, r, n̂)
    K_6 = calc_stiffness_element(H_6, r, n̂)
    return hcat(K_1, K_2, K_3, K_4, K_5, K_6) ./ tol
end


@testset "point stiffness" begin
    tol = 1.0e-08
    K = zeros(Float64,6,6)
    for k = 1:100
        r = randn(SVector{3,Float64})
        n̂ = normalize(randn(SVector{3,Float64}))
        K_point    = calc_point_spatial_stiffness(r, n̂)
        K_point_fd = calc_stiffness_finite_difference(r, n̂, tol)
        @test sum(abs.(K_point_fd - K_point)) < (100 * tol)
        K .+= K_point
    end
    eig_fact = eigen(K)
    @test all(0.0 .< eig_fact.values)
    @test ((abs(det(eig_fact.vectors)) - 1) < 1.0e-10)
end
