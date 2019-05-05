#
# @testset "cholesky_U_deravitive" begin
#     for n = 1:10
#         for k = 1:100
#             A = randn(MMatrix{n,n,Float64,n^2})
#             A = Symmetric(MMatrix(A * A'))  # needs to be PD
#             Ȧ = randn(MMatrix{n,n,Float64,n^2})
#             Ȧ = Symmetric(MMatrix(Ȧ + Ȧ'))  # doesn't need to be PD
#             U, U⁻¹, U̇ = cholesky_U_deravitive(A, Ȧ)
#             Ȧ_formula = U̇' * U + U' * U̇
#             @test sum(abs.(Ȧ_formula - Ȧ)) < (1.0e-10 * n^2)
#         end
#     end
# end
