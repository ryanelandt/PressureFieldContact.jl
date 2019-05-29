
@testset "volume-volume" begin

    ### Box properties
    box_rad = 0.05
    box_density = 400.0
    w_z_0 = 1.14
    i_prop_compliant = InertiaProperties(box_density)
    c_prop_compliant = ContactProperties(Ē=1.0e6)
    eM_box_compliant = as_tet_eMesh(eMesh_box(box_rad))

    mech_scen = MechanismScenario()  # n_quad_rule=2)

    ### Add planes and boxes
    nt_plane  = add_contact!(     mech_scen, "plane", as_tet_eMesh(eMesh_half_plane()),   c_prop=c_prop_compliant)
    nt_body_1 = add_body_contact!(mech_scen, "box_1", eM_box_compliant, i_prop=i_prop_compliant, c_prop=c_prop_compliant)

    ### Friction
    add_friction_regularize!(mech_scen, nt_plane.id,  nt_body_1.id, μd=0.0, χ=0.0, n_quad_rule=2)

    finalize!(mech_scen)
    set_state_spq!(mech_scen, nt_body_1.joint, trans=SVector(0.0, 0.0, 2*box_rad), w=SVector(0.0, 0.0, w_z_0))

    ### Run forward dynamics
    t_final = 5.0e-0
    rr = Radau_for_MechanismScenario(mech_scen)
    rr.step.h_max = 0.05
    data_time, data_state = integrate_scenario_radau(rr, t_final=t_final)

    @test data_state[end, 9] ≈ w_z_0
end
