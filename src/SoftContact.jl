module SoftContact

using StaticArrays
using Rotations: Quat
using RigidBodyDynamics
using GeometryTypes: HomogenousMesh, Face, Point
using Tri_Tet_Intersections
using Binary_BB_Trees
using NumericalTricks
using LinearAlgebra

const frame_tet_cs = CartesianFrame3D("tet_cs")
const frame_tri_cs = CartesianFrame3D("tri_cs")

include("utility.jl")
include("matrix_transform.jl")
include("material.jl")
include("mesh_cache.jl")
include("mesh_inertia.jl")
include("mesh_body_utility.jl")
include("mechanism_scenario.jl")
include("contact_instructions.jl")
include("extensions.jl")
include("contact_algorithms_non_friction.jl")
include("contact_algorithms_friction.jl")

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
    equiv_volume,

    # mesh_body_utility.jl
    bodyFromMesh!,
    newBodyFromInertia,

    # mechanism_scenario
    BristleID,
    BristleFriction,
    ContactInstructions,
    TractionCache,
    TypedQuadTriCache,
    TypedTriTetCache,
    TypedElasticBodyBodyCache,
    TypedMechanismScenario,
    MechanismScenario,

    # contact_instructions.jl
    addContactRigidCompliant!,

    # extensions.jl

    # contact_algorithms_non_friction.jl
    calcXd!,
    refreshJacobians!,
    forceAllElasticIntersections!,
    calcTriTetIntersections!,
    refreshBodyBodyCache!,
    integrateOverContactPatch!,
    intersectionTriangulation!,
    fillTractionCacheForTriangle!,
    fillTractionCacheInnerLoop!,
    calcTangentialVelocity,
    addGeneralizedForcesThirdLaw!,
    addGeneralizedForcesExternal!

    # contact_algorithms_friction.jl
    # friction_regularization,
    # friction_bristle,
    # friction_model


end
