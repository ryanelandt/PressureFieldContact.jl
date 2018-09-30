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

include("mesh_inertia.jl")
include("utility.jl")
include("matrix_transform.jl")
include("material.jl")
include("mesh_cache.jl")
include("mesh_body_utility.jl")
include("mechanism_scenario.jl")
include("contact_instructions.jl")
include("extensions.jl")
include("contact_algorithms_non_friction.jl")
include("contact_algorithms_friction.jl")

export
    frame_tet_cs,
    frame_tri_cs,

    # mesh_inertia.jl
    makeInertiaTensor,
    centroidVolumeCombo,
    equiv_volume,

    # utility.jl
    unPad,
    onePad,
    zeroPad,
    fill_with_nothing!,
    zeroWrench,

    # matrix_transform.jl
    Point4D,
    MatrixTransform,
    getPoint,

    # material.jl
    ContactProperties,
    InertiaProperties,
    calculateExtrensicCompliance,
    calcMutualMu,

    # mesh_cache.jl
    asHomogenousMesh,
    MeshID,
    MeshDict,
    MeshCacheDict,
    MeshCache,
    addBodyMeshCache,

    # mesh_body_utility.jl
    newBodyFromInertia,
    area,
    volume,

    # mechanism_scenario
    BristleID,
    BristleFriction,
    ContactInstructions,
    TractionCache,
    TypedQuadTriCache,
    TypedTriTetCache,
    TypedElasticBodyBodyCache,
    TypedMechanismScenario,
    makePaths,
    MechanismScenario,

    # contact_instructions.jl
    addContactRigidCompliant!,
    # tune_bristle_friction,

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
    addGeneralizedForcesExternal!,

    # contact_algorithms_friction.jl
    regularized_friction,
    find_contact_pressure_center,
    normal_wrench


end
