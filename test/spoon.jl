using CoordinateTransformations: Translation, LinearMap
using StaticArrays
using Rotations: RotZ, RotX, SPQuat
using Binary_BB_Trees
using LinearAlgebra
using RigidBodyDynamics
using MeshCat
using ColorTypes: RGBA, RGB
using MeshCatMechanisms
using SoftContact
using FileIO
# using NumericalTricks
using GeometryTypes: HomogenousMesh


set_zero_subnormals(true)  # optional tweak to make simulation run faster
BLAS.set_num_threads(1)  # optional tweak to make simulation run faster
spoon_path = joinpath(dirname(pathof(SoftContact)), "..", "test", "data", "spoon.obj")

### Create mechanism and temporary structure
my_mechanism = Mechanism(RigidBody{Float64}("world"); gravity=SVector{3,Float64}(0.0, 0.0, -9.8054))  # create empty mechanism
ts = TempContactStruct(my_mechanism)

### Lower box
rad_box = 0.02
lo_name = "box_lo"
lo_h_mesh, lo_tet_mesh = create_volume_box(rad_box, center=SVector(0.0, 0.0, -rad_box))
add_volume_mesh!(ts, root_body(my_mechanism), lo_name, lo_h_mesh, lo_tet_mesh)

### Upper box
up_name = "box_up"
up_joint = Prismatic(SVector{3,Float64}(0.0, 0.0, 1.0))
up_h_mesh, up_tet_mesh = create_volume_box(rad_box)
_, up_joint = add_body_volume_mesh!(ts, up_name, up_h_mesh, up_tet_mesh, InertiaProperties(rho=500.0), joint=up_joint)

### Spoon
mesh_spoon = load(spoon_path)
scale_HomogenousMesh!(mesh_spoon, 0.01)  # make 0.01 of origional size or spoon will be HUGE
# assume that the spoon in a shell of this density and thickness (d) for the purposes of inertia tensor calculation
spoon_i_prop = InertiaProperties(rho=900.0, d=0.0025)
_, spoon_joint = add_body_surface_mesh!(ts, "spoon", mesh_spoon, spoon_i_prop)  # adds body and calculates inertia from the surface mesh

### Friction pairs
add_pair_rigid_compliant_bristle_tune_tri!(ts, "spoon", up_name)
add_pair_rigid_compliant_bristle_tune_tri!(ts, "spoon", lo_name)
mech_scen = MechanismScenario(ts, calcXd!, n_quad_rule=1)  # there will be many intersections; first order is fine

### Set initial condition
# NOTE: All states are zero initialized including deformation state variables
set_state_spq!(mech_scen, spoon_joint, rot=RotX(pi/2), trans=SVector(0.0, 0.0, 0.000001))
set_configuration!(mech_scen.float.state, up_joint, [0.10])
x = get_state(mech_scen)

### Run forward dynamics
h = 0.001
t_final = 1.6
@time data_time, data_state, data_k_iter = integrate_scenario_radau(mech_scen, x*1, t_final=t_final, h_max=h)
println("Finished compiling and running simulation beginning visualization")

### Add meshcat visualizer
if !@isdefined vis
    vis = Visualizer()
    open(vis)
end

### Add meshes to visualizer
mvis = MechanismVisualizer(my_mechanism, vis)
color_pad = RGBA{Float32}(1.0, 0.792, 0.588, 1.0)
lo_h_mesh_32 = HomogenousMesh(vertices=get_h_mesh_vertices_32(lo_h_mesh),
                              faces=get_h_mesh_faces_32(lo_h_mesh),
                              color=color_pad)
setobject!(vis, lo_h_mesh_32)
set_body_mesh_visual!(mvis, mech_scen, "spoon", [0.5, 0.5, 0.5, 1.0])
set_body_mesh_visual!(mvis, mech_scen, up_name, color_pad)

### Move camera
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 27)
settransform!(vis["/Cameras/default"], Translation(0.0, 0.0, 0.04) ∘ LinearMap(RotZ(-π * 0.35)))

### Playback data
sleep(3)  # wait for visualizer to initialize
play_recorded_data(mvis, mech_scen, data_time, data_state, slowdown=5.0)
