
module PressureFieldContact

using DocStringExtensions
using StaticArrays
using Rotations: Quat, Rotation, SPQuat, RotMatrix
using ForwardDiff: Dual
using RigidBodyDynamics
using RigidBodyDynamics.Spatial: vector_to_skew_symmetric, vector_to_skew_symmetric_squared
using GeometryTypes: HomogenousMesh, Face, Point
using ColorTypes: RGBA
using MeshCatMechanisms
using DocStringExtensions
using CoordinateTransformations: Translation

include(joinpath("math_kernel", "MathKernel.jl"))
using .MathKernel

include(joinpath("clip", "Clip.jl"))
using .Clip

include(joinpath("obb", "Binary_BB_Trees.jl"))
using .Binary_BB_Trees

include(joinpath("radau", "Radau.jl"))
using .Radau

include(joinpath("geometry", "Geometry.jl"))
using .Geometry

using LinearAlgebra
using GenericLinearAlgebra


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
    transform!,
    eMesh_half_plane,
    eMesh_sphere,
    eMesh_box,
    eMesh_cylinder,
    as_tet_eMesh,
    as_tri_eMesh,
    n_tri,
    n_tet,
    n_point,
    get_tri,
    get_tet,
    get_point,
    create_swept_mesh,
    f_swept_triv,
    f_swept_circle,

    # math_kernel/MathKernel.jl
    make_pd_gains,
    basic_dh

end
