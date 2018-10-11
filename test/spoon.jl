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
using NumericalTricks
using GeometryTypes: HomogenousMesh


set_zero_subnormals(true)
BLAS.set_num_threads(1)
spoon_path = joinpath(dirname(pathof(SoftContact)), "..", "test", "data", "spoon.obj")

### Box properties
rad_box = 0.02
rad_box_3 = SVector{3,Float64}(1.0, 1.0, 1.0) * rad_box
box_c_prop = ContactProperties(Ē=1.0e6, μ=0.3, d=0.1, χ=1.0)

### Create mechanism and temporary structure
my_mechanism = Mechanism(RigidBody{Float64}("world"); gravity=SVector{3,Float64}(0.0, 0.0, -9.8054))  # create empty mechanism
ts = TempContactStruct(my_mechanism)

### Lower box
lo_name = "box_lo"
lo_mesh, lo_tet_ind, lo_ϵ = outputBoxVolMesh(rad=rad_box_3, center=SVector{3,Float64}(0.0, 0.0, -rad_box))
tet_mesh_lo = TetMesh(lo_mesh.vertices, lo_tet_ind, lo_ϵ, box_c_prop)
add_volume_mesh!(ts, root_body(my_mechanism), lo_name, lo_mesh, tet_mesh_lo)

### Upper box
up_name = "box_up"
up_h_mesh, up_tet_ind, up_ϵ = outputBoxVolMesh(rad=rad_box_3)
up_joint = Prismatic(SVector{3,Float64}(0.0, 0.0, 1.0))
up_tet_mesh = TetMesh(up_h_mesh.vertices, up_tet_ind, up_ϵ, box_c_prop)
add_body_volume_mesh!(ts, up_name, up_h_mesh, up_tet_mesh, InertiaProperties(rho=500.0), up_joint)

### Spoon
mesh_spoon = load(spoon_path)
scale_HomogenousMesh!(mesh_spoon, 0.01)  # make 0.01 of origional size or spoon will be HUGE
joint_spoon = SPQuatFloating{Float64}()
spoon_i_prop = InertiaProperties(rho=900.0, d=0.0025)
add_body_surface_mesh!(ts, "spoon", mesh_spoon, spoon_i_prop, joint_spoon)

### Friction pairs
add_pair_rigid_compliant_bristle_tune_tri!(ts, "spoon", up_name)
add_pair_rigid_compliant_bristle_tune_tri!(ts, "spoon", lo_name)
mech_scen = MechanismScenario(ts, calcXd!, n_quad_rule=1)  # there will be many intersections; first order is fine

### Set initial condition
# NOTE: All states are zero initialized
spoon_rot = components(SPQuat(RotX(pi/2)))
spoon_trans = SVector{3,Float64}(0.0, 0.0, 0.000001)
set_configuration!(mech_scen.float.state, findjoint(my_mechanism, "world_spoon"), vcat(spoon_rot, spoon_trans))  # spoon just off ground
set_configuration!(mech_scen.float.state, findjoint(my_mechanism, "world_" * up_name), [0.10])  # block to 10cm
x = get_state(mech_scen)

### Run forward dynamics
h = 0.001
t_final = 0.2
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
lo_h_mesh_32 = HomogenousMesh(vertices=get_h_mesh_vertices_32(lo_mesh), faces=get_h_mesh_faces_32(lo_mesh), color=color_pad)
setobject!(vis, lo_h_mesh_32)
set_body_mesh_visual!(mvis, mech_scen, "spoon", [0.5, 0.5, 0.5, 1.0])
set_body_mesh_visual!(mvis, mech_scen, up_name, color_pad)

### Move camera
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 27)
settransform!(vis["/Cameras/default"], Translation(0.0, 0.0, 0.04) ∘ LinearMap(RotZ(-π * 0.35)))

### Playback data
sleep(3)  # wait for visualizer to initialize
play_recorded_data(mvis, mech_scen, data_time, data_state, slowdown=5.0)
