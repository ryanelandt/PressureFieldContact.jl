# cd(expanduser("~") * "/.julia/dev/SoftContact/src")

module SoftContact

using StaticArrays
using Rotations: Quat, Rotation, SPQuat
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


const FRAME_ζ¹ = CartesianFrame3D("FRAME_ζ¹")
const FRAME_ζ² = CartesianFrame3D("FRAME_ζ²")
const FRAME_ϕ = CartesianFrame3D("FRAME_ϕ")

include("structs.jl")
include("mesh_inertia.jl")
include("matrix_transform.jl")
include("mesh_cache.jl")
include("mesh_body_utility.jl")
include("temp_structures.jl")
include("mechanism_scenario.jl")
include("extensions.jl")
include("contact_algorithms_non_friction.jl")
include("contact_algorithms_friction.jl")
include("primitive_meshes.jl")
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
    calculateExtrensicCompliance,
    InertiaProperties,
    TetMesh,
    get_tri_mesh,
    get_tet_mesh,
    is_compliant,
    get_Ē,

    # mesh_inertia.jl
    make_volume_mesh_inertia_info,
    make_surface_mesh_inertia_info,
    makeInertiaTensor,
    centroidVolumeCombo,
    equiv_volume,

    # matrix_transform.jl
    Point4D,
    MatrixTransform,
    getPoint,

    # mesh_cache.jl
    asHomogenousMesh,
    MeshID,
    MeshDict,
    MeshCacheDict,
    MeshCache,

    # mesh_body_utility.jl
    newBodyFromInertia,
    outputJointTransform_ParentChild,
    area,
    volume,

    # temp_structures.jl
    BristleID,
    BristleFriction,
    ContactInstructions,
    TempContactStruct,
    addMesh!,
    add_body_volume_mesh!,
    add_volume_mesh!,
    add_body_surface_mesh!,
    add_surface_mesh!,
    add_body_surface!,
    add_body_from_inertia!,
    findmesh,
    findMesh,
    add_pair_rigid_compliant_regularize!,
    add_pair_rigid_compliant!,
    add_pair_rigid_compliant_bristle!,

    # mechanism_scenario
    TractionCache,
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
    verify_bristle_ids!,
    calcXd!,
    refreshJacobians!,
    forceAllElasticIntersections!,
    calcTriTetIntersections!,
    refreshBodyBodyCache!,
    integrate_over_logic!,
    triangle_vertices,
    tetrahedron_vertices_ϵ,
    calc_ζ_transforms,
    find_plane_tet,
    fillTractionCacheForTriangle!,
    fillTractionCacheInnerLoop!,
    calcTangentialVelocity,
    addGeneralizedForcesThirdLaw!,
    addGeneralizedForcesExternal!,

    # contact_algorithms_friction.jl
    regularized_friction,
    find_contact_pressure_center,
    normal_wrench,
    calc_point_spatial_stiffness,
    get_bristle_d0,
    get_bristle_d1,
    # bristle_deformation,
    # bristle_friction!,
    # bristle_friction_no_contact!,

    # primitive_meshes.jl
    create_surface_half_plane,
    create_volume_half_plane,
    create_surface_box,
    create_volume_box,

    # vis_meshcat.jl
    set_body_mesh_visual!,
    HomogenousMesh_32,
    play_recorded_data,

    # example_integrator.jl
    integrate_scenario_radau,

    # utility.jl
    num_partials,
    type_dual,
    fill_with_nothing!,
    mat_mul_SA_bug_circumvent,
    add_h_mesh_color,
    Radau_for_MechanismScenario

end
