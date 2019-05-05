#
# function cholesky_U_deravitive(A::Symmetric{T,MMatrix{N,N,T,N2}}, Ȧ::Symmetric{T,MMatrix{N,N,T,N2}}) where {N,T,N2}
#     # Consider the factorization of a PD matrix A: A = U' * U
#     # The derivative of this factorization is Ȧ = U̇' * U + U' * U̇
#     # This function finds U̇ given Ȧ using the algorithm from: https://homepages.inf.ed.ac.uk/imurray2/pub/16choldiff/choldiff.pdf
#
#     cf = cholesky(A)
#     U = cf.U
#     U⁻¹ = inv(U)
#     term_in = U⁻¹' * Ȧ * U⁻¹
#     @inbounds begin
#         for k = 1:(N+1):N2
#             term_in[k] *= 0.5
#         end
#     end
#     U̇ = UpperTriangular(term_in) * U
#     return U, U⁻¹, U̇
# end
