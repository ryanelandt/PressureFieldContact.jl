using CoordinateTransformations: Translation, LinearMap
using MeshCat
using ColorTypes: RGBA, RGB
using MeshCatMechanisms
using StaticArrays
using LinearAlgebra: BLAS
using RigidBodyDynamics
using PressureFieldContact
using Rotations: RotZ
# using PressureFieldContact.Geometry: eMesh_half_plane, eMesh_box, as_tri_eMesh, as_tet_eMesh
using Test


set_zero_subnormals(true)
BLAS.set_num_threads(1)  # NOTE: comment out this line if using IntelMKL

### Box properties
box_rad = 0.05
box_density = 400.0
c_prop_compliant = ContactProperties(Ē=1.0e6)
i_prop_compliant = InertiaProperties(box_density)
i_prop_rigid     = InertiaProperties(box_density, d=box_rad)
eM_box_rigid     = as_tri_eMesh(eMesh_box(box_rad))
eM_box_compliant = as_tet_eMesh(eMesh_box(box_rad))

mech_scen = MechanismScenario()  # n_quad_rule=2)

### Add planes and boxes
nt_plane  = add_contact!(     mech_scen, "plane", as_tet_eMesh(eMesh_half_plane()),   c_prop=c_prop_compliant)
nt_body_1 = add_body_contact!(mech_scen, "box_1", eM_box_rigid,     i_prop=i_prop_rigid)
nt_body_2 = add_body_contact!(mech_scen, "box_2", eM_box_compliant, i_prop=i_prop_compliant, c_prop=c_prop_compliant)
nt_body_3 = add_body_contact!(mech_scen, "box_3", eM_box_rigid,     i_prop=i_prop_rigid)
nt_body_4 = add_body_contact!(mech_scen, "box_4", eM_box_compliant, i_prop=i_prop_compliant, c_prop=c_prop_compliant)

### Friction
add_friction_regularize!(mech_scen, nt_plane.id,  nt_body_1.id, μd=0.0, χ=2.2, n_quad_rule=2)
add_friction_regularize!(mech_scen, nt_body_1.id, nt_body_2.id, μd=0.2, χ=0.2, n_quad_rule=2)
add_friction_regularize!(mech_scen, nt_body_2.id, nt_body_3.id, μd=0.2, χ=0.2, n_quad_rule=2)
add_friction_regularize!(mech_scen, nt_body_3.id, nt_body_4.id, μd=0.2, χ=0.2, n_quad_rule=2)

finalize!(mech_scen)
set_state_spq!(mech_scen, nt_body_1.joint, trans=SVector(0.0, 0.0,  2*box_rad), w=SVector(0.0, 0.0, 1.0))
set_state_spq!(mech_scen, nt_body_2.joint, trans=SVector(0.0, 0.0,  5*box_rad), w=SVector(0.0, 0.0, 2.0))
set_state_spq!(mech_scen, nt_body_3.joint, trans=SVector(0.0, 0.0,  8*box_rad), w=SVector(0.0, 0.0, 3.0))
set_state_spq!(mech_scen, nt_body_4.joint, trans=SVector(0.0, 0.0, 11*box_rad), w=SVector(0.0, 0.0, 4.0))

if (!@isdefined vis)
    vis = Visualizer();
    (!haskey(ENV, "CI")) && open(vis)
end

### Add meshes to visualizer
color_gray  = RGBA{Float32}(0.5, 0.5, 0.5, 1.0)
color_red   = RGBA{Float32}(1.0, 0.0, 0.0, 1.0)
color_green = RGBA{Float32}(0.0, 1.0, 0.0, 1.0)
color_blue  = RGBA{Float32}(0.0, 0.0, 1.0, 1.0)

mvis = MechanismVisualizer(mech_scen, vis)
set_mesh_visual!(     mvis, mech_scen, nt_plane.id,  color_gray)
set_body_mesh_visual!(mvis, mech_scen, nt_body_1.id, color_red)
set_body_mesh_visual!(mvis, mech_scen, nt_body_2.id, color_blue)
set_body_mesh_visual!(mvis, mech_scen, nt_body_3.id, color_green)
set_body_mesh_visual!(mvis, mech_scen, nt_body_4.id, color_red)
set_configuration!(mech_scen, mvis)

### Run forward dynamics
t_final = 5.0e-0
rr = Radau_for_MechanismScenario(mech_scen)
rr.step.h_max = 0.05
@time data_time, data_state = integrate_scenario_radau(rr, t_final=t_final)
println("Finished compiling and running simulation beginning visualization")

### Move camera
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 7)
settransform!(vis["/Cameras/default"], Translation(0.0, 0.0, 0.30) ∘ LinearMap(RotZ(-π * 0.35)))

### Visualization
if !haskey(ENV, "CI")
    # If experimenting with this script, you can delete all lines below the word "Visualization" except for the one below.
    play_recorded_data(mvis, mech_scen, data_time, data_state, slowdown=1.0)
else
    function run_visualizer()
        play_recorded_data(mvis, mech_scen, data_time, data_state, slowdown=1.0)
        return true
    end

    @testset "visualizer" begin
        @test run_visualizer()
    end
end
