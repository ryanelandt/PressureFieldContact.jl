#
# T = Dual{Nothing,Float64,6}
#
# function partial_length(T::Type)
#     a = T(0.0)
#     return length(a.partials)
# end
#
# rand_partials(N::Int64) = Partials(Tuple(rand(N)))
#
# function rand_sv3_zero(T::Type)
#     (T == Float64) && (return zeros(SVector{3,Float64}))
#     N = partial_length(T)
#     v1 = T(0.0, rand_partials(N))
#     v2 = T(0.0, rand_partials(N))
#     v3 = T(0.0, rand_partials(N))
#     return SVector{3,T}(v1, v2, v3)
# end
#
# function rand_sv3(T::Type)
#     (T == Float64) && (return randn(SVector{3,Float64}))
#     N = partial_length(T)
#     v1 = T(randn(), rand_partials(N))
#     v2 = T(randn(), rand_partials(N))
#     v3 = T(randn(), rand_partials(N))
#     return SVector{3,T}(v1, v2, v3)
# end
#
# is_nan_in_partials(b::Float64) = false
# function is_nan_in_partials(b::Dual{Tag_Type,Float64,N}) where {Tag_Type,N}
#     return any(isnan.(b.partials))
# end
# function is_nan_in_partials(b::SVector{3,T}) where {T}
#     return any(is_nan_in_partials.(b))
# end
#
# @testset "div_by_zero" begin
#     for d_type = [Float64, T]
#         sv3_zero = rand_sv3_zero(d_type)
#         safe_norm_zero = safe_norm(sv3_zero)
#         @test !is_nan_in_partials(safe_norm_zero)
#         @test !is_nan_in_partials(safe_inv_norm(sv3_zero))
#         @test sv3_zero == safe_normalize(sv3_zero)
#         @test !is_nan_in_partials(safe_normalize(sv3_zero))
#         @test !is_nan_in_partials(safe_inv_norm²(sv3_zero))
#         for k = 1:128
#             rv = rand_sv3(d_type)
#             rv_norm = norm(rv)
#             inv_rv_norm = 1 ./ rv_norm
#             rv_normalize = rv .* inv_rv_norm
#             @test norm²(rv) ≈ rv_norm^2
#             @test unsafe_norm(rv) ≈ rv_norm
#             @test safe_norm(rv) ≈ rv_norm
#             @test unsafe_inv_norm(rv) ≈ inv_rv_norm
#             @test safe_inv_norm(rv) ≈ inv_rv_norm
#             @test safe_normalize(rv) ≈ rv_normalize
#             @test unsafe_normalize(rv) ≈ rv_normalize
#             @test safe_inv_norm²(rv) ≈ inv_rv_norm^2
#         end
#         @test zero(T) == safe_scalar_divide(zero(T), rand(T) + 0.5)
#         @test zero(T) == safe_scalar_divide(zero(T), zero(T))
#         @test_throws ErrorException safe_scalar_divide(one(T), zero(T))
#     end
# end
#
# @testset "div by zero type stability test" begin
#     for d_type = [Float64, T]
#         rv = rand_sv3(d_type)
#         rv_0 = rand_sv3_zero(d_type)
#         for fun_name = [norm², unsafe_norm, safe_norm, unsafe_inv_norm, safe_inv_norm, safe_normalize,
#             unsafe_normalize, safe_inv_norm²]
#
#             ans_1 = fun_name(rv)
#             ans_2 = fun_name(rv_0)
#             is_type_stable = typeof(ans_1) == typeof(ans_2)
#             if !is_type_stable
#                 println("for DataType $d_type function $fun_name is not type stable")
#             end
#             @test is_type_stable
#         end
#     end
# end
