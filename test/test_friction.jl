
@testset "bristle friction" begin
    USV = svd(randn(6,6))
    K = Hermitian(Size(6,6)(USV.U * Diagonal(USV.S) * USV.U'))
    sqrt_K, sqrt_K⁻¹ = SoftContact.my_matrix_sqrt_and_inv(K)
    @test (sqrt(K) ≈ sqrt_K)
    @test (inv(sqrt(K).data) ≈ sqrt_K⁻¹)

    A = randn(SMatrix{3,3,Float64,9})
    A_PD = SoftContact.make_3x3_PD(A)
    tol = 1.0e-4
    A_PD == (A + diagm(0=>[1,1,1] * sum(diag(A)) * tol))


end
