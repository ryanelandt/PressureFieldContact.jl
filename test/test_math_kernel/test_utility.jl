
@testset "3_of_6" begin
    a = SVector{6,Float64}(1,2,3,4,5,6)
    @test all(first_3_of_6(a) .== SVector{3,Float64}(1,2,3))
    @test all(last_3_of_6(a) .== SVector{3,Float64}(4,5,6))
end

@testset "utility" begin
    p1 = SVector{3,Float64}(1.0, 2.0, 3.0)
    p2 = SVector{3,Float64}(2.0, 3.0, 4.0)
    @test p2 == weightPoly(p1, p2, 1.0, 0.0)
    @test p1 == weightPoly(p1, p2, 0.0, 1.0)
    @test (p1 + p2) * 0.5 == weightPoly(p1, p2, -0.7, 0.7)

    p, d = make_pd_gains(2*pi, 1.0)
    @test p == 1.0
    @test d == 2.0
    p, d = make_pd_gains(1*pi, 1.0)
    @test p == 4.0
    @test d == 4.0

    for k = 1:8
        pd33 = rand_pd(k)
        @test (isa(pd33, Hermitian{Float64,SMatrix{k,k,Float64,k^2}}))
        @test (0.0 < minimum(eigen(pd33).values))
    end
end

@testset "padding" begin
    v3 = SVector{3,Float64}(7, 8, 9)
    v3_1 = SVector{4,Float64}(v3..., 1)
    v3_0 = SVector{4,Float64}(v3..., 0)

    @test unPad(v3_1) == v3
    @test unPad(SMatrix{1,4,Float64,4}(v3_1...)) == v3
    @test onePad(v3) == v3_1
    @test zeroPad(v3) == v3_0

    M = SMatrix{4,4,Float64,16}(randn(16))
    @test mul_then_un_pad(M, v3_1) == unPad(M * v3_1)
    @test one_pad_then_mul(M, v3) == M * onePad(v3)
end
