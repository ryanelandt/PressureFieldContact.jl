#
# function test_point_bound(p_test, eps_1, clamp_val)
#     n = soft_clamp(p_test - eps_1, clamp_val)
#     z = soft_clamp(p_test,         clamp_val)
#     p = soft_clamp(p_test + eps_1, clamp_val)
#     floating_point_error_allowance = 1.0 + 1.0e-7  # TODO: work on getting this down
#     @test abs(n - z) <= (eps_1 * floating_point_error_allowance)
#     @test abs(p - z) <= (eps_1 * floating_point_error_allowance)
# end
#
# @testset "soft_clamp" begin
#     eps_1 = 1.0e-10
#
#     for scale = (1.0, 2.0)
#         for p_test = [-1.5, -0.5, 0.0, 0.5, 1.5]
#             test_point_bound(p_test * scale, eps_1, 1.0 * scale)
#             test_point_bound(p_test * scale, eps_1, 1.0 * scale)
#         end
#
#         @test -1.0 * scale== soft_clamp(-1.5* scale, 1.0* scale)
#         @test -0.5* scale == soft_clamp(-0.5* scale, 1.0* scale)
#         @test  0.0 * scale== soft_clamp( 0.0* scale, 1.0* scale)
#         @test  0.5 * scale == soft_clamp( 0.5 * scale, 1.0 * scale)
#         @test  1.0 * scale == soft_clamp( 1.5 * scale, 1.0 * scale)
#     end
# end
#
# @testset "fastSigmoid" begin
#     for TYPE_k = (Float64, Dual{Nothing,Float64,12})
#         x_cut⁻¹ = 1.0
#         out = fastSigmoid(TYPE_k( 0.0), x_cut⁻¹)
#         @test zero(TYPE_k) == out
#         @test typeof(out) == TYPE_k
#
#         out = fastSigmoid(TYPE_k( 1.0), x_cut⁻¹)
#         @test  one(TYPE_k) == out
#         @test typeof(out) == TYPE_k
#
#         out = fastSigmoid(TYPE_k(-1.0), x_cut⁻¹)
#         @test -one(TYPE_k) == out
#         @test typeof(out) == TYPE_k
#
#         out = fastSigmoid(TYPE_k( 0.5), x_cut⁻¹)
#         @test zero(TYPE_k) < out < one(TYPE_k)
#         @test typeof(out) == TYPE_k
#
#         out = fastSigmoid(TYPE_k(-0.5), x_cut⁻¹)
#         @test -one(TYPE_k) < out < zero(TYPE_k)
#         @test typeof(out) == TYPE_k
#     end
# end
#
# @testset "poly_approx" begin
#     # fastSoftPlus
#     for center_value = [0.1, 0.25, 0.5, 1.0, 2.0]
#         p = 1.0e-14
#         x = center_value * LinRange(-1.0 - p, -1.0 + p, 1001)
#         y = fastSoftPlus.(x, center_value)
#         @test 0.0 == minimum(y)
#     end
# end
