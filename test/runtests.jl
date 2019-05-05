using Test
using PressureFieldContact
using StaticArrays
using Rotations
using LinearAlgebra
using RigidBodyDynamics
using PressureFieldContact.Tri_Tet_Intersections
using PressureFieldContact.NumericalTricks
using PressureFieldContact.Binary_BB_Trees

set_zero_subnormals(true)
BLAS.set_num_threads(1)  # NOTE: comment out this line if using IntelMKL


include("test_friction.jl")
include("test_exports.jl")
include("test_normal.jl")
include("test_mesh_body_utility.jl")
include("test_vol_vol.jl")
include("boxes.jl")
