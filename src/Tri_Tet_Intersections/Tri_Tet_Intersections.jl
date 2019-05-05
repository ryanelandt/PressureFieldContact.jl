__precompile__(true)

module Tri_Tet_Intersections

using StaticArrays
using ForwardDiff: value
using LinearAlgebra
using NumericalTricks
using RigidBodyDynamics.Spatial


include("poly_eight.jl")
include("utility.jl")
include("quadrature.jl")
include("static_clip.jl")
include("plane_tet_intersection.jl")
include("test_utility.jl")

export
    # poly_eight.jl
    poly_eight,
    zero_small_coordinates,

    # utility.jl
    asMatOnePad,
    asMat,

    # quadrature.jl
    TriTetQuadRule,
    getTriQuadRule,
    getTetQuadRule,

    # static_clip.jl
    poly_eight,
    clip_in_tet_coordinates,

    # plane_tet_intersection.jl
    clip_plane_tet,

    # test_utility.jl
    roll_non_degenerate_tet,
    make_4_sided,
    dist_from_plane,
    project_into_plane

end
