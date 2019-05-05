
@testset "rotations" begin
    qq = Quat(1.0, 0.01, 0.02, 0.03)
    rv = RodriguesVec(qq)
    rv_sv = SVector(rv.sx, rv.sy, rv.sz)
    rv_cheap = cheapRV(qq)
    rel_err = (rv_sv - rv_cheap) ./ rv_sv
    @test all(rel_err .< 1.0e-3)

    t = rand(SVector{4,Float64})
    q = Quat(t..., false)
    @test components(q) == t

    t = rand(SVector{3,Float64})
    s = SPQuat(t...)
    @test components(s) == t

    q = rand(Quat{Float64})
    s = SPQuat(q)
    @test cheapRV(q) ≈ cheapRV(s)

    q_0 = Quat(0.1, 0.2, 0.3, 0.4)
    for k = 1:100
        q_1 = rand(Quat{Float64})
        q_2 = Quat(-components(q_1)..., false)
        q_err_1 = quatErr(q_1, q_0)
        q_err_2 = quatErr(q_2, q_0)
        @test q_err_1 ≈ q_err_1
        q_err_check = components(q_1 * inv(q_0))
        if q_err_check[1] <= 0.0
            @test q_err_check[2:4] ≈ -q_err_1
        else
            @test q_err_check[2:4] ≈ q_err_1
        end
    end
end
