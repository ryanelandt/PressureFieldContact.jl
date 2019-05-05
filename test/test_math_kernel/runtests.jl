# using Test
# using NumericalTricks
# using StaticArrays
# using LinearAlgebra
# using Rotations
# using ForwardDiff: Partials, Dual, value


include("test_exports.jl")
include("test_div_by_zero.jl")
# include("test_matrix_factor_derivatives.jl")
include("test_basic_dh.jl")
include("test_geometry_kernel.jl")
include("test_poly_approx.jl")
include("test_utility.jl")
include("test_vector_projections.jl")
include("test_rotations.jl")
