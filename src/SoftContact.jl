module SoftContact

using RigidBodyDynamics
using StaticArrays
using Tri_Tet_Intersections
using Binary_BB_Trees
using LinearAlgebra
using GeometryTypes: HomogenousMesh, Face, Point


include("material.jl")
include("mesh_cache.jl")
include("mesh_inertia.jl")

export
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
