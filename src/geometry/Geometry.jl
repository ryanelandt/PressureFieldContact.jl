__precompile__(true)

module Geometry


using DocStringExtensions
using GeometryTypes: HomogenousMesh
using Rotations
using LinearAlgebra
using StaticArrays
using DataStructures
using NearestNeighbors
using ..Binary_BB_Trees
using ..NumericalTricks


include("mesh.jl")
include("mesh_create_rot_sym.jl")
include("mesh_create_swept.jl")
include("extensions.jl")
include("blob_types.jl")
include("top_down.jl")
include("util.jl")



export
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
    crop_mesh,
    extrude_mesh,
    output_eMesh_cylinder,

    # mesh_create_rot_sym.jl
    obj_from_point_sequence,

    # mesh_create_swept.jl
    create_swept_mesh,
    f_swept_circle,
    f_swept_triv,

    # blob_types.jl
    blob,
    eMesh_to_tree,

    # util.jl
    calc_min_max,
    calc_obb

    # extensions.jl


end
