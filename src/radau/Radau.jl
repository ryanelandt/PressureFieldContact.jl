module Radau


using ForwardDiff
using StaticArrays
using LinearAlgebra
using DelimitedFiles

# using Polynomials: Poly, polyval, polyder, coeffs
# using PolynomialRoots: roots
# using GenericLinearAlgebra: eigvals
# using GenericSVD: svd

include("load_table_from_file.jl")
include("radau_struct.jl")
include("radau_functions.jl")
include("radau_utilities.jl")
include("radau_solve.jl")
include("interpolate.jl")
include("adaptive.jl")

export
    # load_table_from_file.jl



    # radau_struct.jl
    RadauIntegrator,
    RadauTable,
    RadauStep,
    RadauRule,

    # radau_functions.jl
    calcJacobian!,
    updateInvC!,
    zeroFill!,
    updateFX!,
    calcEw!,
    updateStageX!,

    # radau_utilities.jl
    makeRadauIntegrator,
    get_X_final,
    get_exponent,
    print_exit_flag,

    # radau_solve.jl
    solveRadau,

    # adaptive.jl
    calc_x̂_minus_x,
    calc_x_err_norm,
    update_x_err_norm!,
    update_h!,
    calc_h_new_estimate_1,

    # interpolate.jl
    interpolate_output!

end
