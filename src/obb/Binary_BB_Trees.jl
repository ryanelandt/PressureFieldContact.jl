__precompile__(true)

module Binary_BB_Trees

using Printf
using GeometryTypes: HomogenousMesh
using Rotations
using LinearAlgebra
using StaticArrays
using DataStructures
using Statistics
using NearestNeighbors
using ..NumericalTricks

include("box_types.jl")
include("mesh.jl")
include("util.jl")
include("vector_cache.jl")
include("tree_types.jl")
include("bb_intersection.jl")
include("blob_types.jl")
include("top_down.jl")
include("obb_construction.jl")
include("mesh_create_rot_sym.jl")
include("mesh_create_swept.jl")
include("extensions.jl")


export
    # box_types.jl
    BoundingBox,
    OBB,

    # mesh.jl
    Tri,
    Tet,
    eMesh,
    as_tet_eMesh,
    as_tri_eMesh,
    vertex_pos_for_tri_ind,
    vertex_pos_for_tet_ind,
    n_point,
    n_tri,
    n_tet,
    get_tri,
    get_tet,
    get_point,
    verify_mesh,
    eMesh_transform!,
    invert!,
    mesh_inplace_rekey!,
    mesh_remove_unused_points!,
    delete_triangles!,
    sub_div_mesh,
    mesh_repair!,
    output_eMesh_half_plane,
    output_eMesh_sphere,
    output_eMesh_box,

    # util.jl
    calc_min_max,
    calc_obb,
    crop_mesh,
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

    # blob_types.jl
    blob,
    eMesh_to_tree,

    # top_down.jl

    # obb_construction.jl
    fit_tet_obb,
    fit_tri_obb,

    # mesh_create_rot_sym.jl
    obj_from_point_sequence,

    # mesh_create_swept.jl
    create_swept_mesh,
    f_swept_circle,
    f_swept_triv

    # extensions.jl

end # module
