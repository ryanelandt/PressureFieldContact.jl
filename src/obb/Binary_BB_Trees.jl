__precompile__(true)

module Binary_BB_Trees


using DocStringExtensions
# using Printf
# using GeometryTypes: HomogenousMesh
using Rotations
using LinearAlgebra
using StaticArrays
# using DataStructures
# using Statistics
# using NearestNeighbors
using ..NumericalTricks

include("box_types.jl")
# include("mesh.jl")
include("util.jl")
include("vector_cache.jl")
include("tree_types.jl")
include("bb_intersection.jl")
# include("blob_types.jl")
# include("top_down.jl")
include("obb_construction.jl")
include("extensions.jl")


export
    # box_types.jl
    BoundingBox,
    OBB,

    # util.jl
    calc_min_max,
    calc_obb,
    sort_so_big_ϵ_last,

    # vector_cache.jl
    VectorCache,
    expand!,
    returnNext,
    addCacheItem!,

    # tree_types.jl
    bin_BB_Tree,
    TT_Cache,
    update_TT_Cache!,
    is_leaf,
    is_not_leaf,
    treeDepth,
    leafNumber,
    tree_tree_intersect,
    extractData,

    # bb_intersection.jl
    BB_BB_intersect,

    # obb_construction.jl
    fit_tet_obb,
    fit_tri_obb

    # extensions.jl

end # module
