
@testset "bristle friction" begin
    value_2 = 10.0
    v3_bad = [1.0, 0.0, 0.0]  # zeros on diagonal will fail cholesky decomposition
    K = MMatrix{6,6,Float64,36}(diagm(0=>[v3_bad..., value_2 * v3_bad...]))
    SoftContact.make_stiffness_PD!(K)
    v3_good = 1.0e-4 .+ v3_bad
    K_check = MMatrix{6,6,Float64,36}(diagm(0=>[v3_good..., value_2 * v3_good...]))
    @test all(K .== K_check)
end
