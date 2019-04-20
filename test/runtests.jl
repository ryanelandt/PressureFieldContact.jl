using Test
using SoftContact
using StaticArrays
using Tri_Tet_Intersections
using Rotations
using LinearAlgebra
using NumericalTricks
using RigidBodyDynamics
using Binary_BB_Trees

set_zero_subnormals(true)
BLAS.set_num_threads(1)  # NOTE: comment out this line if using IntelMKL


include("test_exports.jl")
include("test_friction.jl")
include("test_normal.jl")
include("test_mesh_body_utility.jl")
include("test_vol_vol.jl")
