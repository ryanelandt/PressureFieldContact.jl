
using PressureFieldContact: calc_clamped_piecewise, traction

# TODO: remove this test when GenericLinearAlgebra tags a patch that includes PR #53
using GenericLinearAlgebra
@testset "convergence criterion when diagonal has zeros" begin
    M1 = Hermitian(zeros(3, 3))
    M1[1,1] = 0.01
    M2 = Hermitian(zeros(3, 3))
    M2[3, 3] = 0.01
    @test eigvals(M1) == GenericLinearAlgebra._eigvals!(copy(M1))
    @test eigvals(M2) == GenericLinearAlgebra._eigvals!(copy(M2))
    @test eigen(M1).values == GenericLinearAlgebra._eigen!(copy(M1)).values
    @test eigen(M2).values == GenericLinearAlgebra._eigen!(copy(M2)).values
end

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

function verify_fric_calc_bri(ins::Bristle, T̄s::SVector{3,T}, p_dA::T) where {T}
	mag_T̄s = norm(T̄s)
	T_μs = ins.μs * T̄s / mag_T̄s * p_dA
	T_μd = ins.μd * T̄s / mag_T̄s * p_dA
	if mag_T̄s <= ins.μs
		return T̄s * p_dA
	elseif mag_T̄s <= ins.T̄s_μs
		return T_μs
	elseif ins.T̄s_μd <= mag_T̄s
		return T_μd
	else
		w_μd = (mag_T̄s - ins.T̄s_μs) / (ins.T̄s_μd - ins.T̄s_μs)
		w_μs = 1 - w_μd
		return T_μs * w_μs + T_μd * w_μd
	end
end

function verify_fric_calc_reg(ins::Regularized, vel_t::SVector{3,T}, p_dA::T) where {T}
	mag_vel_t = norm(vel_t)
	T_μs = -ins.μs * vel_t / mag_vel_t * p_dA
	T_μd = -ins.μd * vel_t / mag_vel_t * p_dA
	if mag_vel_t <= ins.v_c
		return -ins.μs * vel_t / ins.v_c * p_dA
	elseif mag_vel_t <= ins.v_μs
		return T_μs
	elseif ins.v_μd <= mag_vel_t
		return T_μd
	else
		w_μd = (mag_vel_t - ins.v_μs) / (ins.v_μd - ins.v_μs)
		w_μs = 1 - w_μd
		return T_μs * w_μs + T_μd * w_μd
	end
end

v_c = 1.0e-4
μs_ = 1.1
μd_ = 0.3
ins_bri = Bristle(BristleID(1), τ=0.01, k̄=1000.0, μs=μs_, μd=μd_, magic=1.0e-2)
T̂s_ = normalize(SVector(1.0, 2.0, 3.0))
p_dA_ = 0.133
@testset "bristle" begin
	for v_mag = LinRange(0.0, 4 * μs_, 100)
		T̄s = v_mag * T̂s_
		Tc = traction(ins_bri, T̄s, p_dA_)
		@test Tc ≈ verify_fric_calc_bri(ins_bri, T̄s, p_dA_)
	end
end

ins_reg = Regularized(v_c, μs=μs_, μd=μd_)
v_dir_ = normalize(SVector(1.0, 2.0, 3.0))
p_dA_ = 0.133
@testset "regularized" begin
	for v_mag = LinRange(0.0, 4 * v_c, 100)
		vel_t = v_mag * v_dir_
		Tc = traction(ins_reg, vel_t, p_dA_)
		@test Tc ≈ verify_fric_calc_reg(ins_reg, vel_t, p_dA_)
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

using ForwardDiff: Dual

@testset "decompose and factor K" begin
	for T = (Float64, Dual{Nothing,Float64,9})
		s = spatialStiffness{T}()
		M = Diagonal([1.0, 1.0, 1.0, 1000, 1000, 1000])
		s.K.data .= M * rand_pd(6) * M
		magic = 1.0e-2
		PressureFieldContact.decompose_K!(s, magic)
		K̄ = inv(s.K̄⁻¹_sqrt)^2
		t_1, t_2 = PressureFieldContact.calc_trace_12(K̄)
		@test t_1 ≈ (t_2 * magic^2)
		S = inv(s.S⁻¹)
		@test s.K ≈ S * K̄ * S
	end
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
	S = inv(s.S⁻¹)
	K² = S * inv(s.K̄⁻¹_sqrt * s.K̄⁻¹_sqrt) * S
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
