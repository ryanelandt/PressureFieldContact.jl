# using Test
# using Radau
# using LinearAlgebra

# using Polynomials: Poly, polyval, polyder, coeffs
# using PolynomialRoots: roots
# include("test_T_lambda.jl")


include("test_exports.jl")
include("basic_test.jl")
include("test_robertson.jl")
# include("test_generate_butcher_table.jl")
include("test_time_dep.jl")



# TODO: add tests for calculating of Eigencalculatino of T and λ

# using GenericLinearAlgebra: eigvals
# using GenericSVD: svd
