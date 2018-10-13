
using CoordinateTransformations: Translation, LinearMap
using StaticArrays
using LinearAlgebra: BLAS
using RigidBodyDynamics
using MeshCat
using ColorTypes: RGBA, RGB
using MeshCatMechanisms
using SoftContact
using GeometryTypes: HomogenousMesh
using Binary_BB_Trees: get_h_mesh_faces_32, get_h_mesh_vertices_32
using Rotations: RotZ


set_zero_subnormals(true)
BLAS.set_num_threads(1)

### Box properties
box_rad = 0.05
box_shell = box_rad / 5

### Create mechanism and temporary structure
my_mechanism = Mechanism(RigidBody{Float64}("world"); gravity=SVector{3,Float64}(0.0, 0.0, -9.8054))  # create empty mechanism
ts = TempContactStruct(my_mechanism)

### Half plane
plane_w = 1.0
plane_name = "plane"
plane_h_mesh, plane_tet_mesh = create_volume_half_plane(plane_w, μ=0.0, χ=0.1)
add_volume_mesh!(ts, root_body(my_mechanism), plane_name, plane_h_mesh, plane_tet_mesh)

### Box 1
box_1_name = "box_1"
box_1_h_mesh = create_surface_box(box_rad)
_, box_1_joint = add_body_surface_mesh!(ts, box_1_name, box_1_h_mesh, InertiaProperties(rho=400.0, d=box_shell))

### Box 2
box_2_name = "box_2"
box_2_h_mesh, box_2_tet_mesh = create_volume_box(box_rad, χ=0.03)
_, box_2_joint = add_body_volume_mesh!(ts, box_2_name, box_2_h_mesh, box_2_tet_mesh, InertiaProperties(rho=400.0))

### Box 3
box_3_name = "box_3"
box_3_h_mesh = create_surface_box(box_rad)
_, box_3_joint = add_body_surface_mesh!(ts, box_3_name, box_3_h_mesh, InertiaProperties(rho=400.0, d=box_shell))

### Box 2
box_4_name = "box_4"
box_4_h_mesh, box_4_tet_mesh = create_volume_box(box_rad, χ=0.03)
_, box_4_joint = add_body_volume_mesh!(ts, box_4_name, box_4_h_mesh, box_4_tet_mesh, InertiaProperties(rho=400.0))

### Friction pairs
add_pair_rigid_compliant_regularize!(ts, box_1_name, plane_name)
add_pair_rigid_compliant_regularize!(ts, box_1_name, box_2_name)
add_pair_rigid_compliant_regularize!(ts, box_3_name, box_2_name)
add_pair_rigid_compliant_regularize!(ts, box_3_name, box_4_name)
mech_scen = MechanismScenario(ts, calcXd!, n_quad_rule=2)

### Set initial condition  # NOTE: All states are zero initialized
set_state_spq!(mech_scen, box_1_joint, trans=SVector(0.0, 0.0, 2*box_rad), w=SVector(0.0, 0.0, 1.0))
set_state_spq!(mech_scen, box_2_joint, trans=SVector(0.0, 0.0, 5*box_rad), w=SVector(0.0, 0.0, 2.0))
set_state_spq!(mech_scen, box_3_joint, trans=SVector(0.0, 0.0, 8*box_rad), w=SVector(0.0, 0.0, 3.0))
set_state_spq!(mech_scen, box_4_joint, trans=SVector(0.0, 0.0, 11*box_rad), w=SVector(0.0, 0.0, 4.0))
x = get_state(mech_scen)

### Add meshcat visualizer
if !@isdefined vis
    vis = Visualizer()
    open(vis)
end

### Add meshes to visualizer
color_gray  = RGBA{Float32}(0.5, 0.5, 0.5, 1.0)
color_red   = RGBA{Float32}(1.0, 0.0, 0.0, 1.0)
color_green = RGBA{Float32}(0.0, 1.0, 0.0, 1.0)
color_blue  = RGBA{Float32}(0.0, 0.0, 1.0, 1.0)
mvis = MechanismVisualizer(my_mechanism, vis)
plane_h_mesh_32 = HomogenousMesh_32(plane_h_mesh, color=color_gray)
setobject!(vis[:plane], plane_h_mesh_32)
set_body_mesh_visual!(mvis, mech_scen, box_1_name, color_red)
set_body_mesh_visual!(mvis, mech_scen, box_2_name, color_blue)
set_body_mesh_visual!(mvis, mech_scen, box_3_name, color_green)
set_body_mesh_visual!(mvis, mech_scen, box_4_name, color_red)
set_configuration!(mech_scen, mvis, x)

### Run forward dynamics
h = 0.004
t_final = 5.0e-0
@time data_time, data_state, data_k_iter = integrate_scenario_radau(mech_scen, x*1, t_final=t_final, h_max=h)
println("Finished compiling and running simulation beginning visualization")

### Move camera
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 7)
settransform!(vis["/Cameras/default"], Translation(0.0, 0.0, 0.30) ∘ LinearMap(RotZ(-π * 0.35)))

### Playback data
play_recorded_data(mvis, mech_scen, data_time, data_state, slowdown=2.0)
