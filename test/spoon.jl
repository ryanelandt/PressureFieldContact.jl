
using CoordinateTransformations: Translation, LinearMap
using StaticArrays
using Rotations: RotZ, RotX, SPQuat
using Binary_BB_Trees
using LinearAlgebra  # TODO: is this used?
using RigidBodyDynamics
using MeshCat
using ColorTypes: RGBA, RGB
using MeshCatMechanisms
using PressureFieldContact
using FileIO

set_zero_subnormals(true)  # optional tweak to make simulation run faster
BLAS.set_num_threads(1)  # optional tweak to make simulation run faster (Disable if using MKL)
spoon_path = joinpath(dirname(pathof(PressureFieldContact)), "..", "test", "data", "spoon.obj")

### Create mechanism and temporary structure
my_mechanism = Mechanism(RigidBody{Float64}("world"); gravity=SVector{3,Float64}(0.0, 0.0, -9.8054))  # create empty mechanism
ts = TempContactStruct(my_mechanism)

### Properties + Names
name_lo = "box_lo"
name_up = "box_up"
name_spoon = "spoon"
rad_box = 0.02
i_prop_box = InertiaProperties(rho=500.0)
c_prop_box = ContactProperties(Ē=1.0e6, d=rad_box)

### Lower Box
eM_box_compliant = output_eMesh_box(rad_box, SVector{3,Float64}(0.0, 0.0, -rad_box))
add_contact!(ts, name_lo, eM_box_compliant, c_prop_box)

### Upper Box
joint_up = Prismatic(SVector{3,Float64}(0.0, 0.0, 1.0))
body_up, joint_up = add_body_contact!(ts, name_up, eM_box_compliant, c_prop_box, i_prop_box, joint_type=joint_up)

### Spoon
mesh_spoon = load(spoon_path)
scale_HomogenousMesh!(mesh_spoon, 0.01)  # make 0.01 of origional size or spoon will be HUGE
eMesh_spoon = eMesh(mesh_spoon, nothing, nothing)
# assume that the spoon in a shell of this density and thickness (d) for the purposes of inertia tensor calculation
i_prop_spoon = InertiaProperties(rho=900.0, d=0.0025)
_, spoon_joint = add_body_contact!(ts, "spoon", eMesh_spoon, nothing, i_prop_spoon)  # adds body and calculates inertia from the surface mesh

### Friction
add_friction_bristle!(ts, name_spoon, name_lo, μ=0.2, χ=0.2)
add_friction_bristle!(ts, name_spoon, name_up, μ=0.2, χ=0.2)
mech_scen = MechanismScenario(ts, calcXd!, n_quad_rule=1)  # there will be many intersections; first order is fine

### Set initial condition
# NOTE: All states are zero initialized including deformation state variables
set_state_spq!(mech_scen, spoon_joint, rot=RotX(pi/2), trans=SVector(0.0, 0.0, 0.000001))
set_configuration!(mech_scen.float.state, joint_up, [0.10])
x = get_state(mech_scen)

### Run forward dynamics
t_final = 1.6
rr = Radau_for_MechanismScenario(mech_scen)
data_time, data_state = integrate_scenario_radau(rr, mech_scen, x*1, t_final=t_final)
println("Finished compiling and running simulation beginning visualization")

### Add meshcat visualizer
if !@isdefined vis
    vis = Visualizer()
    open(vis)
end

### Add meshes to visualizer
mvis = MechanismVisualizer(my_mechanism, vis)
color_pad = RGBA{Float32}(1.0, 0.792, 0.588, 1.0)
set_mesh_visual!(mvis, mech_scen, name_lo, color_pad)
set_body_mesh_visual!(mvis, mech_scen, name_spoon, [0.5, 0.5, 0.5, 1.0])
set_body_mesh_visual!(mvis, mech_scen, name_up, color_pad)

### Move camera
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 27)
settransform!(vis["/Cameras/default"], Translation(0.0, 0.0, 0.04) ∘ LinearMap(RotZ(-π * 0.35)))

### Playback data
play_recorded_data(mvis, mech_scen, data_time, data_state, slowdown=1.0)
