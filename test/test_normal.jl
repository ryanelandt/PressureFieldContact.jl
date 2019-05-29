
@testset "normal wrench" begin
    for k_quad_rule = (1, 2)
        p_pos = (0.1, 0.2)

        ### Box properties
        box_rad = 0.05
        i_prop_rigid     = InertiaProperties(400.0, d=0.09)
        c_prop_compliant = ContactProperties(Ē=1.0e9)  # , d=box_rad)

        ### Create mechanism and temporary structure
        mech_scen = MechanismScenario() # n_quad_rule=k_quad_rule)

        name_plane = "plane"
        eM_plane = eMesh_half_plane()
        add_contact!(mech_scen, name_plane, as_tet_eMesh(eM_plane), c_prop=c_prop_compliant)

        name_box = "box"
        eM_box = as_tri_eMesh(eMesh_box(box_rad))
        transform!(eM_box, SVector{3,Float64}(0.0, 0.0, box_rad))
        body_box, joint_box = add_body_contact!(mech_scen, name_box, eM_box, i_prop=i_prop_rigid)

        m_id_plane = find_mesh_id(mech_scen, name_plane)
        m_id_box = find_mesh_id(mech_scen, name_box)
        add_friction_bristle!(mech_scen, m_id_box, m_id_plane, μd=0.3, χ=0.6, k̄=1.0e6, τ=0.03, n_quad_rule=k_quad_rule)
        finalize!(mech_scen)

        pene = 0.1 * box_rad
        set_state_spq!(mech_scen, joint_box, trans=SVector{3,Float64}(p_pos..., -pene))

        x = get_state(mech_scen)
        calcXd(x, mech_scen)

        b = mech_scen.float.bodyBodyCache
        wrench = normal_wrench(b)
        wrench = as_static_vector(wrench)
        check = b.Ē * pene / 1.0 * box_rad^2*4
        f3 = SVector(0.0, 0.0, check)
        a3 = cross(SVector(p_pos..., 0.0), f3)
        check = vcat(a3,f3)
        @test wrench ≈ check

        # wrench, p_center = normal_wrench_patch_center(b)
        # wrench = as_static_vector(wrench)
        # check = b.Ē * pene / 1.0 * box_rad^2*4
        # @test check ≈ wrench[6]
        # @test all((p_pos .* 0.999) .< p_center.v[1:2] .< (p_pos .* 1.001))
    end
end
