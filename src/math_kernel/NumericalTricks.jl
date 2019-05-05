__precompile__(true)

module NumericalTricks

using StaticArrays
using LinearAlgebra
using Rotations
using ForwardDiff: Dual, value


include("geometry_kernel.jl")
include("poly_approx.jl")
include("div_by_zero.jl")
include("rotations.jl")
include("vector_projections.jl")
include("utility.jl")
include("basic_dh.jl")
# include("matrix_factor_derivatives.jl")


export
  # geometry_kernel.jl
  area,
  centroid,
  triangleNormal,
  triangle_area,
  vector_area,
  volume,

  # # poly_approx.jl
  # fastSoftPlus,
  # fastSigmoid,
  # soft_clamp,
  # # smooth_c1_ramp,
  # # smooth_c2_ramp,

  # # div_by_zero.jl
  # safe_normalize,
  # unsafe_normalize,
  # safe_inv_norm²,
  # safe_inv_norm,
  # safe_norm,
  # norm²,
  # safe_scalar_divide,
  # unsafe_norm,
  # unsafe_inv_norm,

  # rotations.jl
  quatErr,
  cheapRV,
  components,

  # vector_projections.jl
  vec_proj,
  vec_sub_vec_proj,
  a_dot_one_pad_b,

  # utility.jl
  first_3_of_6,
  last_3_of_6,
  getTop,
  unPad,
  onePad,
  zeroPad,
  mul_then_un_pad,
  one_pad_then_mul,
  weightPoly,
  make_pd_gains,
  rand_pd,

  # basic_dh.jl
  basic_dh,
  dh_R_t,
  povray_12,
  dh_vector_mul

  # # matrix_factor_derivatives.jl
  # cholesky_U_deravitive

end
