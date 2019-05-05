# function fastSoftPlus(x::T, k::Float64=1.0) where {T}
#     # This function is designed to ALWAYS produce non-negative output
#     (x <= -k) && (return zero(T))
#     (k <= x) && (return x)
#     return (k + x)^3 * (3*k - x) / (16*k^3)
# end
#
# function fastSigmoid(x_::T, x_sat⁻¹::Float64=6.0) where {T}
#     # using SymPy
#     # using LinearAlgebra
#     #
#     # @syms x a0 a1 a2 a3 a4 a5 x_cut real=true
#     #
#     # p_0 = a1*x + a3*x^3 + a5*x^5
#     # p_1 = diff(p_0, x, 1)
#     # p_2 = diff(p_1, x, 1)
#     #
#     # e1 = subs(p_0 - 1, x, x_cut)
#     # e2 = subs(p_1 - 0, x, x_cut)
#     # e3 = subs(p_2 - 0, x, x_cut)
#     #
#     # A,B = SymPy.sympy["linear_eq_to_matrix"]([e1, e2, e3], a1, a3, a5)
#     # c = A \ B
#     #
#     # p = dot(c, [x^k for k = 1:2:5])
#     # println("    poly = ", p)
#
#     x = x_ * x_sat⁻¹
#     x = clamp(x, -one(T), one(T))
#     x2 = x * x
#     term = muladd(0.375, x2, -1.25)
#     return x * muladd(term, x2, +1.875)
# end
#
# function soft_clamp(x::T, bound::T) where {T}
#     # using LinearAlgebra
#     # using SymPy
#     #
#     # @syms x a0 a1 a2 a3 a4 a5 real=true
#     #
#     # p_0 = a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5
#     # p_1 = diff(p_0, x, 1)
#     # p_2 = diff(p_1, x, 1)
#     #
#     # e1 = subs(p_0 - 1//2, x, 1//2)
#     # e2 = subs(p_1 - 1, x, 1//2)
#     # e3 = subs(p_2 - 0, x, 1//2)
#     #
#     # e4 = subs(p_0 - 1, x, 3//2)
#     # e5 = subs(p_1 - 0, x, 3//2)
#     # e6 = subs(p_2 - 0, x, 3//2)
#     #
#     # A,B = SymPy.sympy["linear_eq_to_matrix"]([e1, e2, e3, e4, e5, e6], a0, a1, a2, a3, a4, a5)
#     # c = A \ B
#     #
#     # p = dot(c, [x^k for k = 0:5])
#     # println("    poly = ", p)
#     x̄ = safe_scalar_divide(x, bound)
#     return soft_clamp_normalized(x̄, x, bound)
# end
#
# function soft_clamp_normalized(x̄::T, x::T, bound::T) where {T}
#     (x̄ < 0) && return -soft_clamp_normalized(-x̄, -x, bound)
#     (1.5 <= x̄) && (return bound)
#     (x̄ <= 0.5) && (return x)
#     x2 = x̄ * x̄
#     poly = x2 * (x2 * 0.5 - 2 * x̄ + Float64(9/4)) + Float64(5/32)
#     return bound * poly
# end
#
# #############################################################
#
# # smooth_c1_ramp(v1::T, v2::T) where {T} = v2 * smooth_c1_ramp(safe_scalar_divide(v1, v2))
# # function smooth_c1_ramp(t::T) where {T}
# #     # y(t) = 1.5 * t - 0.5 * t^3 -- satisfies the following derivative conditions
# #     #
# #     #  0 = y(0) = ÿ(0) = ẏ(1) = ẏ(-1)
# #     #  1 = y(1)
# #     # -1 = y(-1)
# #
# #     if 1 <= abs(t)
# #         return sign(t) * one(T)
# #     else
# #         return 0.5 * t * (3 - t * t)
# #     end
# # end
# #
# # smooth_c2_ramp(v1::T, v2::T) where {T} = v2 * smooth_c2_ramp(safe_scalar_divide(v1, v2))
# # function smooth_c2_ramp(t::T) where {T}
# #     # y(t) = 1.875 * t - 1.25 * t^3 + 0.375 * t^5 -- satisfies the following derivative conditions
# #     #
# #     #  0 = y(0) = ÿ(0) = y⁴(0) = ẏ(1) = ẏ(-1) = ÿ(1) = ÿ(-1)
# #     #  1 = y(1)
# #     # -1 = y(-1)
# #
# #     if 1 <= abs(t)
# #         return sign(t) * one(T)
# #     else
# #         t² = t * t
# #         return 0.125 * t * (15 + t² * (3 * t² - 10))
# #     end
# # end
#
#
# # function fastSigmoid(x::T) where {T}
# #     # This function assumes that x is non-negative.
# #     (0.16666666666666666 <= x) && (return one(T))
# #     x2 = x * x
# #     x3 = x * x2
# #     x5 = x3 * x2
# #     return 11.25 * x - 270.0 * x3 + 2916.0 * x5
# # end
#
#
# # # y = d + c 1 x + b 1 x^2 + a 1 x^3  -- 0
# # # y = 0 + c 1   + b 2 x   + a 3 x^2  -- 1
# # # y = 0 +   0   + b 2     + a 6 x    -- 2
# #
# # using UnicodePlots
# # using LinearAlgebra
# #
# #
# # function the_coe(poly_len::Int64, n_to_go::Int64)
# #     if n_to_go == 0
# #         return 1
# #     else
# #         return (poly_len - 1) * the_coe(poly_len - 1, n_to_go - 1)
# #     end
# # end
# # function the_formula(t::Float64, d::Int64, k::Int64)
# #     the_deg = k - 1 - d
# #     if the_deg == 0
# #         t_pow = 1.0
# #     elseif the_deg <= -1
# #         t_pow = 0.0
# #     else
# #         t_pow = t ^ the_deg
# #     end
# #     the_coe(k, d) * t_pow
# # end
# # expand_poly_at(c, d::Int64, t::Float64) = [c[k] * the_formula(t, d, k) for k = 1:n]
# # eval_poly_at(c, d::Int64, t::Float64) = sum(expand_poly_at(c, d, t))
# # function add_thing(c_, d::Int64, t, b_)
# #     push!(A, expand_poly_at(c_, d, t))
# #     push!(b, b_)
# #     return nothing
# # end
# #
# # n = 6
# # c0 = [1 for k = 1:n]
# # A = Vector{Vector{Float64}}()
# # b = zeros(0)
# # t_1 = 1.0
# # add_thing(c0, 0, 0.0, 0.0)
# # add_thing(c0, 2, 0.0, 0.0)
# # add_thing(c0, 0, t_1, sign(t_1))
# # add_thing(c0, 1, t_1, 0.0)
# # add_thing(c0, 4, 0.0, 0.0)  # second order only
# # add_thing(c0, 2, t_1, 0.0)  # second order only
# #
# # A = vcat(A'...)
# # c = A \ b
# #
# # t_space = collect(LinRange(0.0, t_1, 100))
# # y_eval = eval_poly_at.([c], 0, t_space)
# #
# # myPlot = lineplot(t_space, y_eval)
# # println(myPlot)
# # println(c)
# #
# # t_space = collect(LinRange(-1.5, 1.5, 100))
# # myPlot = lineplot(t_space, smooth_c2_ramp.(t_space))
# # println(myPlot)
