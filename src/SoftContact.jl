
module SoftContact

using Printf
using StaticArrays
using Rotations: Quat, Rotation, SPQuat, RotMatrix
using ForwardDiff: Dual
using RigidBodyDynamics
using RigidBodyDynamics.Spatial: vector_to_skew_symmetric
using GeometryTypes: HomogenousMesh, Face, Point
using ColorTypes: RGBA
using MeshCatMechanisms
using Tri_Tet_Intersections
using CoordinateTransformations: Translation
using Binary_BB_Trees
using NumericalTricks
using LinearAlgebra
using Radau
using GenericLinearAlgebra


const FRAME_ζ¹ = CartesianFrame3D("FRAME_ζ¹")
const FRAME_ζ² = CartesianFrame3D("FRAME_ζ²")
const FRAME_ϕ = CartesianFrame3D("FRAME_ϕ")

include("structs.jl")
include("mesh_inertia.jl")
include("matrix_transform.jl")
include("mesh_body_utility.jl")
include("temp_structures.jl")
include("mechanism_scenario.jl")
include("extensions.jl")
include("contact_algorithms_non_friction.jl")
include("contact_algorithms_friction.jl")
include("contact_algorithms_normal.jl")
include("vis_meshcat.jl")
include("example_integrator.jl")
include("utility.jl")

export
    FRAME_ζ¹,
    FRAME_ζ²,
    FRAME_ϕ,

    # structs.jl
    MeshInertiaInfo,
    ContactProperties,
    eTree,
    InertiaProperties,
    get_tree_tet,
    get_tree_tri,
    get_c_prop,
    get_ind_tri,
    get_ind_tet,
    get_ϵ,
    is_compliant,
    is_rigid,
    get_Ē,
    MeshID,
    MeshDict,
    MeshCacheDict,
    MeshCache,

    # mesh_inertia.jl
    make_volume_mesh_inertia_info,
    make_surface_mesh_inertia_info,
    makeInertiaInfo,
    makeInertiaTensor,
    centroidVolumeCombo,
    equiv_volume,

    # matrix_transform.jl
    Point4D,
    MatrixTransform,
    getPoint,

    # mesh_body_utility.jl
    newBodyFromInertia,
    outputJointTransform_ParentChild,
    # area,
    # volume,

    # temp_structures.jl
    BristleID,
    Bristle,
    Regularized,
    ContactInstructions,
    TempContactStruct,
    addMesh!,
    add_body_contact!,
    make_eTree_obb,
    add_contact!,
    add_body!,
    add_body_from_inertia!,  # this needs to be here
    add_friction_regularize!,
    add_friction!,
    add_friction_bristle!,

    # mechanism_scenario
    TractionCache,
    calc_p_dA,
    spatialStiffness,
    TypedElasticBodyBodyCache,
    TypedMechanismScenario,
    makePaths,
    MechanismScenario,
    num_partials,
    num_x,
    type_dual,
    get_state,
    set_state_spq!,

    # extensions.jl
    principal_value!,

    # contact_algorithms_non_friction.jl
    calcXd!,
    calcXd,
    refreshJacobians!,
    normal_wrench,
    forceAllElasticIntersections!,
    calcTriTetIntersections!,
    refreshBodyBodyCache!,
    # integrate_over_logic!,
    addGeneralizedForcesThirdLaw!,

    # contact_algorithms_friction.jl
    set_s_from_Δʷ,

    # contact_algorithms_normal.jl
    normal_wrench,

    # vis_meshcat.jl
    set_body_mesh_visual!,
    set_mesh_visual!,
    HomogenousMesh_32,
    asHomogenousMesh,
    play_recorded_data,

    # example_integrator.jl
    integrate_scenario_radau,

    # utility.jl
    num_partials,
    type_dual,
    fill_with_nothing!,
    # mat_mul_SA_bug_circumvent,
    # add_h_mesh_color,
    Radau_for_MechanismScenario,
    get_bristle_d0,
    get_bristle_d1,
    as_static_vector,
    find_mesh,
    find_mesh_id

end
