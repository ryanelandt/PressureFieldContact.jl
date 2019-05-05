
function create_box_and_plane(n_quad_rule::Int64, tang_force_coe::Float64=0.0, v_tol::Float64=NaN)
    box_rad = 0.05
    i_prop_compliant = InertiaProperties(400.0)
    i_prop_rigid     = InertiaProperties(400.0, d=0.09)
    Ē = 1.0e9
    c_prop_compliant = ContactProperties(Ē=Ē)
    eM_box_rigid     = as_tri_eMesh(output_eMesh_box(box_rad))
    eM_box_compliant = output_eMesh_box(box_rad)

    ### Create mechanism and temporary structure
	mag_grav = 9.8054
    mech_scen = MechanismScenario(n_quad_rule=n_quad_rule)

	eM_plane = output_eMesh_half_plane(1.0)
    nt_plane = add_contact!(mech_scen, "plane", as_tet_eMesh(eM_plane), c_prop=c_prop_compliant)

    eM_box_1 = output_eMesh_box(box_rad * ones(SVector{3,Float64}))
	eMesh_transform!(eM_box_1, SVector(0, 0, box_rad))
	nt_box_1 = add_body_contact!(mech_scen, "box_1", as_tri_eMesh(eM_box_1), i_prop=i_prop_rigid)

	μ = 0.3
	is_regularize = !isnan(v_tol)
	if is_regularize
		add_friction_regularize!(mech_scen, nt_plane.id, nt_box_1.id, μ=μ, v_tol=v_tol)
		vel_box = SVector(0.0, v_tol, 0.0)
	else
		add_friction_bristle!(mech_scen, nt_plane.id, nt_box_1.id, μ=μ, k̄=1.0e4, τ=0.03)
		vel_box = SVector(0.0, 0.0, 0.0)
	end
	finalize!(mech_scen)
	mass_time_grav = mag_grav * nt_box_1.body.inertia.mass
	area_box_face =  4 * box_rad^2
	pene = mass_time_grav / (Ē * area_box_face)
	set_state_spq!(mech_scen, nt_box_1.joint, vel=vel_box, trans=SVector(0.0, 0.0, -pene))
	mech_scen.τ_ext[5] = mass_time_grav * μ * tang_force_coe
	return mech_scen
end


i_box_y_vel = 11
v_tol = 1.0e-4

@testset "spatial spring friction" begin
	for n_quad_rule = 1:2
		# Test that velocity of box is positive when pushed with 1.5 the friction strength
		mech_scen = create_box_and_plane(n_quad_rule, 1.5)
		rr = Radau_for_MechanismScenario(mech_scen)
		data_time, data_state = integrate_scenario_radau(rr, t_final=0.1)
		@test 0.00001 < data_state[end, i_box_y_vel]

		# Test that velocity of the box is zero when pushed with 0.5 the friction strength
		mech_scen = create_box_and_plane(n_quad_rule, 0.5)
		rr = Radau_for_MechanismScenario(mech_scen)
		data_time, data_state = integrate_scenario_radau(rr, t_final=1.0)
		@test abs(data_state[end, i_box_y_vel]) < 1.0e-12
	end
end

@testset "regularized friction" begin
	for n_quad_rule = 1:2
		# Test that the box slows down when traveling at v_tol and a force less than friction strength is applied
		mech_scen = create_box_and_plane(n_quad_rule, 0.999, v_tol)
		@test calcXd(get_state(mech_scen), mech_scen)[i_box_y_vel] < 0.0

		# Test that the box speeds up when traveling at v_tol and a force less than friction strength is applied
		mech_scen = create_box_and_plane(n_quad_rule, 1.001, v_tol)
		@test 0.0 < calcXd(get_state(mech_scen), mech_scen)[i_box_y_vel]
	end
end


# @testset "transform_stiffness" begin
#     f2 = CartesianFrame3D()
#     f3 = CartesianFrame3D()
#     I3 = @SMatrix([9.0 2.0 1.0; 2.0 8.0 3.0; 1.0 3.0 7.0])
#     inertia_pre = SpatialInertia{Float64}(f2, I3, randn(SVector{3,Float64}), 1.5)
#     xform = Transform3D(f2, f3, rand(RotMatrix{3,Float64}), rand(SVector{3,Float64}))
#     i66 = SMatrix(transform(inertia_pre, xform))
#
#     s = spatialStiffness{Float64}()
#     s.K.data .= SMatrix(inertia_pre)
#     PressureFieldContact.transform_stiffness!(s, xform)
#
#     @test i66 ≈ s.K
# end

@testset "calc_K_sqrt⁻¹" begin
    s = spatialStiffness{Float64}()
	K_orig = Matrix(rand_pd(6))
    s.K.data .= K_orig
    PressureFieldContact.calc_K_sqrt⁻¹!(s)
    @test s.K⁻¹_sqrt ≈ sqrt(inv(K_orig))
end

@testset "spatial spring friction" begin
    ### Box properties
    box_rad = 0.05
    i_prop_compliant = InertiaProperties(400.0)
    i_prop_rigid     = InertiaProperties(400.0, d=0.09)
    Ē = 1.0e9
    c_prop_compliant = ContactProperties(Ē=Ē)
    eM_box_rigid     = as_tri_eMesh(output_eMesh_box(box_rad))
    eM_box_compliant = output_eMesh_box(box_rad)

    ### Create mechanism and temporary structure
    mech_scen = MechanismScenario(n_quad_rule=2)

    name_part = "part"
	eM_box = output_eMesh_half_plane(1.0)
    nt_part = add_body_contact!(mech_scen, name_part, as_tri_eMesh(eM_box), i_prop=i_prop_rigid)

    name_hol_1 = "hol_1"
    hol_rad = 0.2 * box_rad
    eM_hol_1 = output_eMesh_box(hol_rad * ones(SVector{3,Float64}))
	eMesh_transform!(eM_hol_1, SVector{3,Float64}(0.0, 0.0, hol_rad))  # box_rad + hol_rad, 0.0, 0.0))

    joint_1 = Prismatic(SVector{3,Float64}(0.0, 0.0, 1.0))
    nt_hol_1 = add_body_contact!(mech_scen, name_hol_1, as_tet_eMesh(eM_hol_1), c_prop=c_prop_compliant, i_prop=i_prop_compliant, joint=joint_1)

    τ = 0.03
    k̄ = 1.0e6
    add_friction_bristle!(mech_scen, nt_part.id, nt_hol_1.id, μ=0.3, χ=0.6, k̄=k̄, τ=τ)

    finalize!(mech_scen)  # , calcXd!, n_quad_rule=2)
	mech_scen = mech_scen

    pene = hol_rad * 0.001
    set_configuration!(mech_scen.float.state, nt_hol_1.joint, [-pene])
    x = get_state(mech_scen)

	# #######################
	# if !@isdefined vis  ### Add meshcat visualizer
	# 	vis = Visualizer()
	# 	open(vis)
	# end
	# color_gray  = RGBA{Float32}(0.5, 0.5, 0.5, 1.0)
	# color_blue  = RGBA{Float32}(0.0, 0.0, 1.0, 1.0)
	# mvis = MechanismVisualizer(my_mechanism, vis)
	# set_mesh_visual!(mvis, mech_scen, nt_hol_1.id, color_blue)
	# set_mesh_visual!(mvis, mech_scen, nt_part.id,  color_gray)
	# set_configuration!(mech_scen, mvis, x)
	# #######################

    ### Test 1 -- see if stiffness is what is expected
    calcXd(x, mech_scen)
    K_ana = (hol_rad^2) * 4 * k̄ * (Ē * (pene / hol_rad))
	s = mech_scen.float.bodyBodyCache.spatialStiffness
	K² = inv(s.K⁻¹_sqrt * s.K⁻¹_sqrt)
	K_44 = K²[4,4]
	K_55 = K²[5,5]
	@test K_44 ≈ K_55
    @test 0.99 * K_ana < K_55 < 1.01 * K_ana

    # ### Test 2 -- see if stiffness was decomposed correctly
    # tm = mech_scen.float
    # b  = tm.bodyBodyCache
    # spaStiff = b.spatialStiffness
	#
    # ### Test 3 -- see if spring force calculated both ways agrees
    # c_ins = mech_scen.ContactInstructions[1]
    # Δ² = (rand(SVector{6,Float64}) + 0.5) * 1.0e-8
	#
    # s = PressureFieldContact.set_s_from_Δʷ(mech_scen, c_ins, Δ²)
    # mech_scen.float.s.segments[c_ins.FrictionModel.BristleID] .= s
    # _, wrench²_fric = PressureFieldContact.bristle_wrench_in_world(tm,  c_ins)
    # wrench²_fric = as_static_vector(wrench²_fric)
	# wrench²_fric_2 = - K² * Δ²
    # @test (0.999999 <  dot(normalize(wrench²_fric), normalize(wrench²_fric_2)))
end
