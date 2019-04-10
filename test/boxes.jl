using CoordinateTransformations: Translation, LinearMap
using StaticArrays
using LinearAlgebra: BLAS
using RigidBodyDynamics
using MeshCat
using ColorTypes: RGBA, RGB
using MeshCatMechanisms
using SoftContact
using Rotations: RotZ
using Binary_BB_Trees: output_eMesh_half_plane, output_eMesh_box, as_tri_eMesh, scale!


set_zero_subnormals(true)
BLAS.set_num_threads(1)  # NOTE: comment out this line if using IntelMKL

### Box properties
box_rad = 0.05
i_prop_compliant = InertiaProperties(400.0)
i_prop_rigid = InertiaProperties(400.0, d=0.09)
c_prop_compliant = ContactProperties(Ē=1.0e6, d=box_rad)
eM_box_rigid = as_tri_eMesh(output_eMesh_box(box_rad))
eM_box_compliant = output_eMesh_box(box_rad)

### Create mechanism and temporary structure
my_mechanism = Mechanism(RigidBody{Float64}("world"); gravity=SVector{3,Float64}(0.0, 0.0, -9.8054))  # create empty mechanism
ts = TempContactStruct(my_mechanism)

### Names
name_plane = "plane"
name_box_1 = "box_1"
name_box_2 = "box_2"
name_box_3 = "box_3"
name_box_4 = "box_4"

### Plane
d_plane = 1.0
eM_half = output_eMesh_half_plane(d_plane)
m_id_plane = add_contact!(ts, name_plane, eM_half, ContactProperties(Ē=1.0e6, d=d_plane))

### Boxes
nt_body_1 = add_body_contact!(ts, name_box_1, eM_box_rigid, nothing, i_prop_rigid)
nt_body_2 = add_body_contact!(ts, name_box_2, eM_box_compliant, c_prop_compliant, i_prop_compliant)
nt_body_3 = add_body_contact!(ts, name_box_3, eM_box_rigid, nothing, i_prop_rigid)
nt_body_4 = add_body_contact!(ts, name_box_4, eM_box_compliant, c_prop_compliant, i_prop_compliant)

### Friction
add_pair_rigid_compliant_regularize!(ts, m_id_plane, nt_body_1.mesh_id, μ=0.0, χ=2.2)
add_pair_rigid_compliant_regularize!(ts, nt_body_1.mesh_id, nt_body_2.mesh_id, μ=0.2, χ=0.2)
add_pair_rigid_compliant_regularize!(ts, nt_body_2.mesh_id, nt_body_3.mesh_id, μ=0.2, χ=0.2)
add_pair_rigid_compliant_regularize!(ts, nt_body_3.mesh_id, nt_body_4.mesh_id, μ=0.2, χ=0.2)

mech_scen = MechanismScenario(ts, calcXd!, n_quad_rule=2)
set_state_spq!(mech_scen, nt_body_1.joint, trans=SVector(0.0, 0.0, 2*box_rad), w=SVector(0.0, 0.0, 1.0))
set_state_spq!(mech_scen, nt_body_2.joint, trans=SVector(0.0, 0.0, 5*box_rad), w=SVector(0.0, 0.0, 2.0))
set_state_spq!(mech_scen, nt_body_3.joint, trans=SVector(0.0, 0.0, 8*box_rad), w=SVector(0.0, 0.0, 3.0))
set_state_spq!(mech_scen, nt_body_4.joint, trans=SVector(0.0, 0.0, 11*box_rad), w=SVector(0.0, 0.0, 4.0))
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
set_mesh_visual!(mvis, mech_scen, m_id_plane, color_gray)
set_body_mesh_visual!(mvis, mech_scen, nt_body_1.mesh_id, color_red)
set_body_mesh_visual!(mvis, mech_scen, nt_body_2.mesh_id, color_blue)
set_body_mesh_visual!(mvis, mech_scen, nt_body_3.mesh_id, color_green)
set_body_mesh_visual!(mvis, mech_scen, nt_body_4.mesh_id, color_red)
set_configuration!(mech_scen, mvis, x)

### Run forward dynamics
t_final = 5.0e-0
rr = Radau_for_MechanismScenario(mech_scen)
data_time, data_state = integrate_scenario_radau(rr, mech_scen, x*1, t_final=t_final)
println("Finished compiling and running simulation beginning visualization")

### Move camera
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 7)
settransform!(vis["/Cameras/default"], Translation(0.0, 0.0, 0.30) ∘ LinearMap(RotZ(-π * 0.35)))

### Playback data
play_recorded_data(mvis, mech_scen, data_time, data_state, slowdown=1.0)
