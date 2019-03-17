
@testset "bristle friction" begin
    ### Box properties
    box_rad = 0.05
    i_prop_compliant = InertiaProperties(400.0)
    i_prop_rigid     = InertiaProperties(400.0, d=0.09)
    c_prop_compliant = ContactProperties(Ē=1.0e9, d=box_rad)
    eM_box_rigid     = as_tri_eMesh(output_eMesh_box(box_rad))
    eM_box_compliant = output_eMesh_box(box_rad)

    ### Create mechanism and temporary structure
    my_mechanism = Mechanism(RigidBody{Float64}("world"); gravity=SVector{3,Float64}(0.0, 0.0, -9.8054))  # create empty mechanism
    ts = TempContactStruct(my_mechanism)

    name_part = "part"
    eM_box = output_eMesh_box(box_rad .* SVector{3,Float64}(1.0, 1.0, 0.25))
    dh_transform_mesh!(eM_box, basic_dh(RotZ(pi/4)))
    body_part, joint_part = add_body_contact!(ts, name_part, eM_box, c_prop_compliant, i_prop_compliant)

    rad_inner = 0.010
    rad_outer = 0.020
    r1 = (rad_outer + rad_inner) / 2
    r2 = (rad_outer - rad_inner) / 2
    name_hol_1 = "hol_1"

    lr = LinRange(pi/4, 3*pi/4, 5)
    eM_hol_1 = create_swept_mesh(θ->f_swept_circle(r1, θ), lr, r2)

    x3 = 1.15 * box_rad

    dh_transform_mesh!(eM_hol_1, basic_dh(SVector{3,Float64}(0.0, x3, 0.0)))
    add_contact!(ts, name_hol_1, eM_hol_1, c_prop_compliant)

    m_id_part = find_mesh_id(ts, name_part)
    m_id_hol_1 = find_mesh_id(ts, name_hol_1)
    add_pair_rigid_compliant_bristle!(ts, m_id_part, m_id_hol_1, μ=0.3, χ=0.6, k̄=1.0e6, τ=0.03)

    mech_scen = MechanismScenario(ts, calcXd!, n_quad_rule=2)
    x = get_state(mech_scen)

    c_ins = mech_scen.ContactInstructions[1]
    tm = mech_scen.float
    b = tm.bodyBodyCache



    SoftContact.force_single_elastic_intersection!(mech_scen, tm, c_ins)  # populate the cache

    ### Test 1
    # see if stiffness was decomposed correctly
    s = b.spatialStiffness
    SoftContact.decompose_stiffness!(s)
    K_back = s.S * (s.Ū' * s.Ū) * s.S
    @test K_back ≈ b.spatialStiffness.K

    ### Test 2
    # see if spring force calculated both ways agrees
    δ = (rand(SVector{6,Float64}) + 0.5) * 1.0e-8
    Δ = SoftContact.calc_Δ(s, δ)
    mech_scen.float.s.segments[BristleID(1)] .= Δ
    _, wrenchʷ_fric = SoftContact.bristle_wrench_in_world(tm, c_ins)
    wrenchʷ_fric = as_static_vector(wrenchʷ_fric)
    wrenchʷ_fric_2 = - s.Kʷ * δ
    @test all(wrenchʷ_fric ≈ wrenchʷ_fric_2)
end
