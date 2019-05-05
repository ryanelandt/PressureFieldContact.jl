
module PressureFieldContact

# using Printf
using StaticArrays
using Rotations: Quat, Rotation, SPQuat, RotMatrix
using ForwardDiff: Dual
using RigidBodyDynamics
using RigidBodyDynamics.Spatial: vector_to_skew_symmetric, vector_to_skew_symmetric_squared
using GeometryTypes: HomogenousMesh, Face, Point
using ColorTypes: RGBA
using MeshCatMechanisms
using CoordinateTransformations: Translation

include(joinpath("math_kernel", "NumericalTricks.jl"))
using .NumericalTricks

include(joinpath("Tri_Tet_Intersections", "Tri_Tet_Intersections.jl"))
using .Tri_Tet_Intersections

include(joinpath("obb", "Binary_BB_Trees.jl"))
using .Binary_BB_Trees

include(joinpath("radau", "Radau.jl"))
using .Radau

using LinearAlgebra
using GenericLinearAlgebra
using DocStringExtensions

# const FRAME_ζ¹ = CartesianFrame3D("FRAME_ζ¹")
# const FRAME_ζ² = CartesianFrame3D("FRAME_ζ²")
# const FRAME_ϕ = CartesianFrame3D("FRAME_ϕ")

include("structs.jl")
include("body_inertia.jl")
# include("matrix_transform.jl")
include("mechanism_scenario.jl")
include("extensions.jl")
include("contact_algorithms_non_friction.jl")
include("contact_algorithms_friction.jl")
include("contact_algorithms_normal.jl")
include("vis_meshcat.jl")
include("example_integrator.jl")
include("utility.jl")

export
    # FRAME_ζ¹,
    # FRAME_ζ²,
    # FRAME_ϕ,

    # structs.jl
    MeshInertiaInfo,
    ContactProperties,
    # eTree,
    InertiaProperties,
    get_c_prop,
    get_ind_tri,
    get_ind_tet,
    get_ϵ,
    get_Ē,
    MeshID,
    MeshDict,
    MeshCacheDict,
    MeshCache,

    # mesh_inertia.jl
    makeInertiaInfo,
    makeInertiaTensor,
    centroidVolumeCombo,
    equiv_volume,

    # # matrix_transform.jl
    # Point4D,
    # MatrixTransform,
    # getPoint,

    # mesh_body_utility.jl
    newBodyFromInertia,
    outputJointTransform_ParentChild,

    # mechanism_scenario
    BristleID,
    Bristle,
    Regularized,
    ContactInstructions,
    TractionCache,
    calc_p_dA,
    spatialStiffness,
    TypedElasticBodyBodyCache,
    TypedMechanismScenario,
    DiscreteControl,
    MechanismScenario,
    finalize!,
    num_partials,
    num_x,
    type_dual,
    get_state,
    set_state_spq!,
    addMesh!,
    add_body_contact!,
    # make_eTree_obb,
    add_contact!,
    add_body!,
    add_body_from_inertia!,  # this needs to be here
    add_friction_regularize!,
    add_friction_bristle!,

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
    addGeneralizedForcesThirdLaw!,

    # contact_algorithms_friction.jl

    # contact_algorithms_normal.jl
    normal_wrench,

    # vis_meshcat.jl
    set_body_mesh_visual!,
    set_mesh_visual!,
    HomogenousMesh_32,
    play_recorded_data,

    # example_integrator.jl
    integrate_scenario_radau,

    # utility.jl
    num_partials,
    type_dual,
    fill_with_nothing!,
    Radau_for_MechanismScenario,
    get_bristle_d0,
    get_bristle_d1,
    as_static_vector,
    find_mesh,
    find_mesh_id,

    # obb/Binary_BB_Trees.jl
    output_eMesh_half_plane,
    output_eMesh_sphere,
    output_eMesh_box,
    as_tet_eMesh,
    as_tri_eMesh,
    n_tri,
    n_tet,
    n_point,
    get_tri,
    get_tet,
    get_point,

    # math_kernel/NumericalTricks.jl
    make_pd_gains

end
