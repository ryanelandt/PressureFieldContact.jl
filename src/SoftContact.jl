module SoftContact

using RigidBodyDynamics
using StaticArrays
using Tri_Tet_Intersections
using Binary_BB_Trees
using LinearAlgebra
using GeometryTypes: HomogenousMesh, Face, Point


const frame_tet_cs = CartesianFrame3D("tet_cs")
const frame_tri_cs = CartesianFrame3D("tri_cs")

include("utility.jl")
include("matrix_transform.jl")
include("material.jl")
include("mesh_cache.jl")
include("mesh_inertia.jl")

export
    frame_tet_cs,
    frame_tri_cs,

    # utility.jl
    unPad,
    onePad,
    zeroPad,
    fill_with_nothing!,

    # matrix_transform.jl
    Point4D,
    MatrixTransform,
    getPoint,

    # material.jl
    ContactMaterialProperties,
    InertiaMaterialProperties,
    MaterialProperties,
    calculateExtrensicCompliance,
    calcMutualMu,

    # mesh_cache.jl
    RawMeshCache,
    asHomogenousMesh,
    MeshID,
    MeshDict,
    MeshCacheDict,
    MeshCache,

    # mesh_inertia.jl
    makeInertiaTensor,
    centroidVolumeCombo,
    equiv_volume

end
