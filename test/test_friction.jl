
@testset "bristle friction" begin
    ### Box properties
    box_rad = 0.05
    i_prop_compliant = InertiaProperties(400.0)
    i_prop_rigid     = InertiaProperties(400.0, d=0.09)
    Ē = 1.0e9
    c_prop_compliant = ContactProperties(Ē=Ē, d=box_rad)
    eM_box_rigid     = as_tri_eMesh(output_eMesh_box(box_rad))
    eM_box_compliant = output_eMesh_box(box_rad)

    ### Create mechanism and temporary structure
    my_mechanism = Mechanism(RigidBody{Float64}("world"); gravity=SVector{3,Float64}(0.0, 0.0, -9.8054))  # create empty mechanism
    ts = TempContactStruct(my_mechanism)

    name_part = "part"
    eM_box = output_eMesh_box(box_rad .* SVector{3,Float64}(1.0, 1.0, 20.0))
    dh_transform_mesh!(eM_box, basic_dh(box_rad * SVector{3,Float64}(0.0, 0.0, -19.0)))
    body_part, joint_part = add_body_contact!(ts, name_part, as_tri_eMesh(eM_box), nothing, i_prop_rigid)

    name_hol_1 = "hol_1"
    hol_rad = 0.2 * box_rad
    eM_hol_1 = output_eMesh_box(hol_rad * ones(SVector{3,Float64}))
    dh_transform_mesh!(eM_hol_1, basic_dh(SVector{3,Float64}(box_rad + hol_rad, 0.0, 0.0)))
    joint_1 = Prismatic(SVector{3,Float64}(-1.0, 0.0, 0.0))
    _, joint_1 = add_body_contact!(ts, name_hol_1, eM_hol_1, c_prop_compliant, i_prop_compliant, joint_type=joint_1)

    name_hol_2 = "hol_2"
    eM_hol_2 = deepcopy(eM_hol_1)
    dh_transform_mesh!(eM_hol_2, basic_dh(RotZ(1.0 * pi)))
    joint_2 = Prismatic(SVector{3,Float64}(+1.0, 0.0, 0.0))
    _, joint_2 = add_body_contact!(ts, name_hol_2, eM_hol_2, c_prop_compliant, i_prop_compliant, joint_type=joint_2)

    m_id_part = find_mesh_id(ts, name_part)
    m_id_hol_1 = find_mesh_id(ts, name_hol_1)
    m_id_hol_2 = find_mesh_id(ts, name_hol_2)
    τ = 0.03
    k̄ = 1.0e6
    add_pair_rigid_compliant_bristle!(ts, m_id_part, m_id_hol_1, μ=0.3, χ=0.6, k̄=k̄, τ=τ)
    add_pair_rigid_compliant_bristle!(ts, m_id_part, m_id_hol_2, μ=0.3, χ=0.6, k̄=k̄, τ=τ)

    mech_scen = MechanismScenario(ts, calcXd!, n_quad_rule=2)
    pene = hol_rad * 0.001
    set_configuration!(mech_scen.float.state, joint_1, [pene])
    set_configuration!(mech_scen.float.state, joint_2, [pene])
    x = get_state(mech_scen)

    ### Test 1 -- see if stiffness is what is expected
    calcXd!(get_state(mech_scen), x, mech_scen)
    K = (hol_rad^2) * 4 * k̄ * (Ē * (pene / hol_rad))
    Kʷ = mech_scen.float.bodyBodyCache.spatialStiffness.Kʷ
    @test 0.99 * K < Kʷ[6,6] < 1.01 * K

    ### Test 2 -- see if stiffness was decomposed correctly
    tm = mech_scen.float
    b = tm.bodyBodyCache
    s = b.spatialStiffness
    SoftContact.decompose_stiffness!(s)
    K_back = s.S * (s.Ū' * s.Ū) * s.S
    @test K_back ≈ b.spatialStiffness.K

    ### Test 3 -- see if spring force calculated both ways agrees
    c_ins = mech_scen.ContactInstructions[1]
    SoftContact.force_single_elastic_intersection!(mech_scen, tm, c_ins)  # populate the cache
    SoftContact.decompose_stiffness!(s)
    δʷ = (rand(SVector{6,Float64}) + 0.5) * 1.0e-8
    Δ = set_Δ_from_δʷ(mech_scen, c_ins, δʷ)
    mech_scen.float.s.segments[c_ins.FrictionModel.BristleID] .= Δ
    _, wrenchʷ_fric = SoftContact.bristle_wrench_in_world(tm,  c_ins)
    wrenchʷ_fric = as_static_vector(wrenchʷ_fric)
    wrenchʷ_fric_2 = - s.Kʷ * δʷ
    @test (0.999999 <  dot(normalize(wrenchʷ_fric), normalize(wrenchʷ_fric_2)))
end
