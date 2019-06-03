
using PressureFieldContact: calc_clamped_piecewise, calc_μ, traction

@testset "clac_clamped_piecewise" begin
	x_1 = 0.3
	x_2 = 0.5
	y_1 = 1.1
	y_2 = 0.1
	@test y_1   ≈ calc_clamped_piecewise(x_1 - 0.1,        x_1, x_2, y_1, y_2)
	@test y_1   ≈ calc_clamped_piecewise(x_1,              x_1, x_2, y_1, y_2)
	@test y_1   ≈ calc_clamped_piecewise(x_1 + 10 * eps(), x_1, x_2, y_1, y_2)
	x_avg = (x_1 + x_2) / 2
	y_avg = (y_1 + y_2) / 2
	@test y_avg ≈ calc_clamped_piecewise(x_avg,            x_1, x_2, y_1, y_2)
	@test y_2   ≈ calc_clamped_piecewise(x_2 - 10 * eps(), x_1, x_2, y_1, y_2)
	@test y_2   ≈ calc_clamped_piecewise(x_2,              x_1, x_2, y_1, y_2)
	@test y_2   ≈ calc_clamped_piecewise(x_2 + 0.1,        x_1, x_2, y_1, y_2)
end

ins_Bristle = Bristle(BristleID(1), τ=0.01, k̄=1.0e4, μs=1.0, μd=0.3)
ins_Reg = Regularized(0.1, μs=1.0, μd=0.1)

@testset "FrictionModel" begin
	@test ins_Bristle.μs < ins_Bristle.T̄s_μ_is_μd
	@test ins_Bristle.δ_max < ins_Bristle.δ_μ_is_0

	@test ins_Reg.v_tol < ins_Reg.v_μ_is_μd
end

@testset "calc_μ" begin
	@test ins_Bristle.μs == calc_μ(ins_Bristle, 0.0, 0.0)
	@test ins_Bristle.μd == calc_μ(ins_Bristle, 0.0, 10.0)
	@test            0.0 == calc_μ(ins_Bristle, 10.0, 0.0)

	@test ins_Reg.μs == calc_μ(ins_Reg, 0.0)
	μ_avg = (ins_Reg.μs + ins_Reg.μd) / 2
	v_avg = (ins_Reg.v_tol + ins_Reg.v_μ_is_μd) / 2
	@test μ_avg ≈ calc_μ(ins_Reg, v_avg)

	@test ins_Reg.μd == calc_μ(ins_Reg, 10.0)
end

@testset "traction" begin
	sv_0 = zeros(SVector{3,Float64})

	for p_dA_ = (0.1, 1.0, 10.0)
		@test sv_0 == traction(ins_Bristle, ins_Bristle.δ_μ_is_0, SVector(0.01, 0.01, 0.01), p_dA_)

		r3 = SVector(0.01, 0.2, 0.3)
		@test r3 * p_dA_ == traction(ins_Bristle, ins_Bristle.δ_max, r3, p_dA_)

		@test normalize(r3) * ins_Bristle.μd * p_dA_ ≈ traction(ins_Bristle, ins_Bristle.δ_max, r3 * 1000, p_dA_)
	end
end

function create_box_and_plane(n_quad_rule::Int64, tang_force_coe::Float64=0.0, v_tol::Float64=NaN)
    box_rad = 0.05
    i_prop_compliant = InertiaProperties(400.0)
    i_prop_rigid     = InertiaProperties(400.0, d=0.09)
    Ē = 1.0e9
    c_prop_compliant = ContactProperties(Ē=Ē)
    eM_box_rigid     = as_tri_eMesh(eMesh_box(box_rad))
    eM_box_compliant = eMesh_box(box_rad)

    ### Create mechanism and temporary structure
	mag_grav = 9.8054
    mech_scen = MechanismScenario()  # n_quad_rule=n_quad_rule)

	eM_plane = eMesh_half_plane(1.0)
    nt_plane = add_contact!(mech_scen, "plane", as_tet_eMesh(eM_plane), c_prop=c_prop_compliant)

    eM_box_1 = eMesh_box(box_rad * ones(SVector{3,Float64}))
	transform!(eM_box_1, SVector(0, 0, box_rad))
	nt_box_1 = add_body_contact!(mech_scen, "box_1", as_tri_eMesh(eM_box_1), i_prop=i_prop_rigid)

	μd = 0.3
	is_regularize = !isnan(v_tol)
	if is_regularize
		add_friction_regularize!(mech_scen, nt_plane.id, nt_box_1.id, μd=μd, v_tol=v_tol, n_quad_rule=n_quad_rule)
		vel_box = SVector(0.0, v_tol, 0.0)
	else
		add_friction_bristle!(mech_scen, nt_plane.id, nt_box_1.id, μd=μd, k̄=1.0e4, τ=0.03, n_quad_rule=n_quad_rule)
		vel_box = SVector(0.0, 0.0, 0.0)
	end
	finalize!(mech_scen)
	mass_time_grav = mag_grav * nt_box_1.body.inertia.mass
	area_box_face =  4 * box_rad^2
	pene = mass_time_grav / (Ē * area_box_face)
	set_state_spq!(mech_scen, nt_box_1.joint, vel=vel_box, trans=SVector(0.0, 0.0, -pene))
	mech_scen.τ_ext[5] = mass_time_grav * μd * tang_force_coe
	return mech_scen
end

i_box_y_vel = 11
v_tol = 1.0e-4

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
    eM_box_rigid     = as_tri_eMesh(eMesh_box(box_rad))
    eM_box_compliant = eMesh_box(box_rad)

    ### Create mechanism and temporary structure
    mech_scen = MechanismScenario()

    name_part = "part"
	eM_box = eMesh_half_plane(1.0)
    nt_part = add_body_contact!(mech_scen, name_part, as_tri_eMesh(eM_box), i_prop=i_prop_rigid)

    name_hol_1 = "hol_1"
    hol_rad = 0.2 * box_rad
    eM_hol_1 = eMesh_box(hol_rad * ones(SVector{3,Float64}))
	transform!(eM_hol_1, SVector{3,Float64}(0.0, 0.0, hol_rad))  # box_rad + hol_rad, 0.0, 0.0))

    joint_1 = Prismatic(SVector{3,Float64}(0.0, 0.0, 1.0))
    nt_hol_1 = add_body_contact!(mech_scen, name_hol_1, as_tet_eMesh(eM_hol_1), c_prop=c_prop_compliant, i_prop=i_prop_compliant, joint=joint_1)

    τ = 0.03
    k̄ = 1.0e6
    add_friction_bristle!(mech_scen, nt_part.id, nt_hol_1.id, μd=0.3, χ=0.6, k̄=k̄, τ=τ, n_quad_rule=2)

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
end

function calc_it(t::SVector{3,Float64})
    box_rad = 0.05
    box_density = 400.0
    c_prop_compliant = ContactProperties(Ē=1.0e6)
    i_prop_rigid     = InertiaProperties(box_density, d=box_rad)
    eM_box_rigid     = as_tri_eMesh(eMesh_box(box_rad))
    mech_scen = MechanismScenario()  # n_quad_rule=2)
    nt_plane  = add_contact!(     mech_scen, "plane", as_tet_eMesh(eMesh_half_plane()),   c_prop=c_prop_compliant)
    nt_body_1 = add_body_contact!(mech_scen, "box_1", eM_box_rigid,     i_prop=i_prop_rigid)
    add_friction_bristle!(mech_scen, nt_plane.id,  nt_body_1.id, μd=1.0, χ=2.2, n_quad_rule=2)
    finalize!(mech_scen)
    set_state_spq!(mech_scen, nt_body_1.joint, trans=t + SVector(0.0, 0.0,  0.99 * box_rad), w=SVector(0.4, 0.3, 1.0))
    calcXd(get_state(mech_scen), mech_scen)
    tm = mech_scen.float
    cop, _ = PressureFieldContact.normal_wrench_cop(tm.bodyBodyCache)
    BF = mech_scen.ContactInstructions[1].FrictionModel
    PressureFieldContact.calc_patch_spatial_stiffness!(mech_scen.float, BF, cop.v)
    return mech_scen.float.bodyBodyCache.spatialStiffness.K, cop.v
end

t = SVector(0.35, 0.10, 0.0)
K_t, cop_t = calc_it(t)
K_0, cop_0 = calc_it(SVector(0.0, 0.0, 0.0))

@testset "spatial stiffness" begin
    @test K_0 ≈ K_t
    @test cop_t ≈ (cop_0 + t)
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
