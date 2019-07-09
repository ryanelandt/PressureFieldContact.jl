
using CoordinateTransformations: Translation, LinearMap
using StaticArrays
using LinearAlgebra: BLAS
using RigidBodyDynamics
using RigidBodyDynamics.Spatial: vector_to_skew_symmetric
using MeshCat
using ColorTypes: RGBA, RGB
using MeshCatMechanisms
using Rotations
using PressureFieldContact
using PressureFieldContact.Radau
using LinearAlgebra


set_zero_subnormals(true)
BLAS.set_num_threads(1)  # NOTE: comment out this line if using IntelMKL


is_bristle = true  # set to false to use REGURALIZED friction between the pencil and the pads

# TUNABLE BRISTLE FRICTION PARAMATERS
k̄ = 8.0e4
τ = 0.01
magic = 1.0e-2

# TUNABLE REGURALIZED FRICTION PARAMATERS
v_tol = 1.0e-5

# Geometric Constants
pad_rad = 0.0035
penci_length = 0.16
penci_rad = 0.0035

mutable struct GripperIns
	z::Float64  # height of the arm
	slip_allowed::Float64  # see comments below GripperIns
	θ_arm::Float64  # desired angle of arm
	Δt::Float64  # when should it look for the next instruction
	t0::Float64  # is set to the current time the first time the instruction is used
	rot_0::Float64  # is set to the current ***pencil*** rotation angle the first time the instruction is used
	function GripperIns(; z::Float64=z_high, slip_allowed::Float64=0.0, θ_arm::Float64=0.0, Δt::Float64=2.0)
		# This object contains information about how the gripper should behave.
		# Much of this behavior is encoded in the paramater called slip_allowed as follows:
		# -Inf    --> no grip (fingers are separated )
		#  0.0    --> no slip allowed (fingers grip tightly)
		#  finite --> finite slip allowed (fingers allow slip until a certain slip is reached)
		#  Inf    --> infinite slip allowed (the fingers grip loosely and allow the pencil to rotate)
		return new(z, slip_allowed, θ_arm, Δt, Inf, Inf)
	end
end

# setpoints for the arm height and rotation anngle
z_high = 1.2 * penci_length
z_low = penci_rad
θ_arm_1 =  pi / 2
θ_arm_2 = -pi / 2

v_gripper_ins = [
	### Pick up Pencil
	# lower gripper to pencil
	GripperIns(Δt=1.2, slip_allowed=-Inf, z=z_low,),
	# grip pencil while lowered
	GripperIns(Δt=0.25, slip_allowed=+Inf, z=z_low),
	# let pencil slip as gripper lifts up
	GripperIns(Δt=1.3, slip_allowed=+Inf),

	### Rotate Arm 90 degrees CCW
	# grip pencil hard to prevent future slip
	GripperIns(Δt=0.5, slip_allowed=0.0),
	# rotate pencil ccw 90 degrees
	GripperIns(Δt=2.0, slip_allowed=0.0, θ_arm=θ_arm_1),
	# loosen grip and let pencil slip
	GripperIns(Δt=1.5, slip_allowed=+Inf, θ_arm=θ_arm_1),

	### Rotate Arm 180 degrees CW
	# grip pencil hard to prevent future slip
	GripperIns(Δt=0.5, slip_allowed=0.0, θ_arm=θ_arm_1),
	# rotate pencil cw 180 degrees
	GripperIns(Δt=2.5, slip_allowed=0.0, θ_arm=θ_arm_2),

	### Control Slipping
	# let pencil slip a bit then stop
	GripperIns(Δt=2.5, slip_allowed=0.2, θ_arm=θ_arm_2),
	# let pencil slip a bit then stop
	GripperIns(Δt=2.5, slip_allowed=0.2, θ_arm=θ_arm_2),
	# let pencil rotate until it falls out of grasp
	GripperIns(Δt=Inf, slip_allowed=Inf, θ_arm=θ_arm_2),
]

function grip_control!(x::Vector{Float64}, m::MechanismScenario, t::Float64)
	joint_to_id(joint::Joint) = mech_scen.float.state.v[joint].indices[1]

	function calc_q̈(j::Joint, q_des::Float64, q̇_des=0.0; q̈_max::Float64=Inf)
		# Calculates the desired acceleration for a PD control law
		q = mech_scen.float.state.q[j][1]
		q̇ = mech_scen.float.state.v[j][1]
		q_err = q - q_des
		q̇_err = q̇ - 0.0
		q̈_des = - gain_d * q̇_err - gain_p * q_err
		q̈_des = clamp(q̈_des, -q̈_max, +q̈_max)
		return q̈_des, q_err, q̇_err
	end

	function get_grip_ins(t::Float64, rot_0::Float64)
		# get the correct item in v_gripper_ins based on how much time has passed
		grip_ins = v_gripper_ins[1]
		(grip_ins.t0 == Inf) && (grip_ins.t0 = t)
		if grip_ins.Δt < (t - grip_ins.t0)  # instruction expired
			popfirst!(v_gripper_ins)
			grip_ins = v_gripper_ins[1]
			grip_ins.t0 = t
			grip_ins.rot_0 = rot_0
		end
		return grip_ins
	end

	gain_p, gain_d = make_pd_gains(0.75, 1.0)  # 0.75s settling time, 1.0 damping ratio

	# the joint is in the named tuple associated with each joint
	j_tra_z = nt_tra_z.joint
	j_rev_y = nt_rev_y.joint
	j_pad_n = nt_pad_n.joint
	j_pad_p = nt_pad_p.joint
	j_penci = nt_penci.joint

	# Each joint has an ID
	id_tra_z = joint_to_id(j_tra_z)
	id_rev_y = joint_to_id(j_rev_y)
	id_pad_n = joint_to_id(j_pad_n)
	id_pad_p = joint_to_id(j_pad_p)

	id_actuate = vcat(id_tra_z, id_rev_y, id_pad_n, id_pad_p)  # only the DoF of the gripper are actuated

	# update dynamics quantaties
	tm = mech_scen.float
	state = tm.state
    copyto!(tm, x)
    H = tm.result.massmatrix
    mass_matrix!(H, state)
    dynamics_bias!(tm.result, state)
    configuration_derivative!(tm.result.q̇, state)
	H = mech_scen.float.result.massmatrix
	H4 = H[id_actuate, id_actuate]
	C4 = mech_scen.float.result.dynamicsbias[id_actuate]

	angle_penci = AngleAxis(rotation(transform_to_root(state, nt_penci.body))).theta
	grip_ins = get_grip_ins(t, angle_penci)

	q̈_tra_z, q_err_tra_z 	= calc_q̈(j_tra_z, grip_ins.z, q̈_max=1.0)
	q̈_rev_y, q_rev_y		= calc_q̈(j_rev_y, grip_ins.θ_arm)
	q̈_pad_n, _ 			= calc_q̈(j_pad_n, 0.0)
	q̈_pad_p, _ 			= calc_q̈(j_pad_p, 0.0)

	q̈ = [q̈_tra_z,	q̈_rev_y, q̈_pad_n, q̈_pad_p]
	m.τ_ext[id_actuate] .= H4 * q̈ + C4  # See Featherstone

	# use the slip_allowed field on the GripperIns to figure out how hard to pinch
	τ_soft = 0.3
	τ_hard = 2.5
	slip_allowed = grip_ins.slip_allowed
	if slip_allowed == -Inf
		τ_grip = 0.0
	elseif slip_allowed == 0.0
		τ_grip = τ_hard
	elseif slip_allowed == Inf
		τ_grip = τ_soft
	elseif isfinite(slip_allowed)
		Δ_ang_penci = abs(angle_penci - grip_ins.rot_0)
		τ_grip = ifelse(Δ_ang_penci < grip_ins.slip_allowed, τ_soft, τ_hard)
	end
	m.τ_ext[id_pad_n] .+= τ_grip  # just add grip forces onto the contant-ignorant forces
	m.τ_ext[id_pad_p] .+= τ_grip
end

### Create mechanism and temporary structure
mech_scen = MechanismScenario(discrete_controller=DiscreteControl(grip_control!, 0.01))
my_mechanism = mech_scen.float.state.mechanism

c_prop_compliant = ContactProperties(Ē=1.0e6)
i_prop_compliant = InertiaProperties(400.0)
i_prop_pad       = InertiaProperties(16000.0)
i_prop_rigid     = InertiaProperties(400.0, d=penci_rad)

eM_plane = eMesh_half_plane()
transform!(eM_plane, 0.6)

eM_pad_n = as_tet_eMesh(eMesh_sphere(pad_rad, 4))
transform!(eM_pad_n, SVector(0.0, -(pad_rad + penci_rad), 0.0))
transform!(eM_pad_n, SMatrix{3,3,Float64,9}(2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0))

name_rev_y = "rev_y"
eM_pad_p = as_tet_eMesh(eMesh_sphere(pad_rad, 4))
transform!(eM_pad_p, SVector(0.0, +(pad_rad + penci_rad), 0.0))
transform!(eM_pad_p, SMatrix{3,3,Float64,9}(2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0))
eM_rev_y = as_tri_eMesh(eMesh_box(pad_rad * SVector(1, 7, 1), pad_rad * SVector(0, -4, 8)))

tip_turn = 0.013
eM_penci = as_tri_eMesh(create_swept_mesh(f_swept_triv, [0.0, tip_turn, penci_length], [0.0, penci_rad, penci_rad], 12, true, rot_half=true))
transform!(eM_penci, SVector(0.0, 0.0, penci_rad))

joint_tra_z = Prismatic(SVector(0.0,  0.0, 1.0))
joint_rev_y = Revolute( SVector(0.0,  1.0, 0.0))
joint_n     = Prismatic(SVector(0.0, +1.0, 0.0))
joint_p     = Prismatic(SVector(0.0, -1.0, 0.0))

nt_plane = add_contact!(mech_scen, "plane", as_tet_eMesh(eM_plane), c_prop=c_prop_compliant)

mii_base = MeshInertiaInfo(1 .* ones(SMatrix{3,3,Float64,9}), zeros(SVector{3,Float64}), 1.0, NaN)
nt_tra_z = add_body_from_inertia!(my_mechanism, "tra_z", mii_base, joint=joint_tra_z)
nt_rev_y = add_body_from_inertia!(my_mechanism, name_rev_y, mii_base, joint=joint_rev_y, body=nt_tra_z.body)
nt_rev_y_id = add_contact!(mech_scen, name_rev_y, eM_rev_y, body=nt_rev_y.body)

nt_pad_n = add_body_contact!(mech_scen, "pad_n", eM_pad_n, c_prop=c_prop_compliant, i_prop=i_prop_pad, joint=joint_n, body=nt_rev_y.body)
nt_pad_p = add_body_contact!(mech_scen, "pad_p", eM_pad_p, c_prop=c_prop_compliant, i_prop=i_prop_pad, joint=joint_p, body=nt_rev_y.body)
nt_penci = add_body_contact!(mech_scen, "name", eM_penci, i_prop=i_prop_rigid)

if is_bristle
	c_ins_pad_n = add_friction_bristle!(mech_scen, nt_penci.id, nt_pad_n.id, μd=0.5, χ=0.6, k̄=k̄, magic=magic, τ=τ)
	c_ins_pad_p = add_friction_bristle!(mech_scen, nt_penci.id, nt_pad_p.id, μd=0.5, χ=0.6, k̄=k̄, magic=magic, τ=τ)
else
	c_ins_pad_n = add_friction_regularize!(mech_scen, nt_penci.id, nt_pad_n.id, μd=0.5, χ=0.6, v_tol=v_tol)
	c_ins_pad_p = add_friction_regularize!(mech_scen, nt_penci.id, nt_pad_p.id, μd=0.5, χ=0.6, v_tol=v_tol)
end
add_friction_regularize!(mech_scen, nt_penci.id, nt_plane.id, μd=0.5, χ=0.6, v_tol=v_tol)
add_friction_regularize!(mech_scen, nt_pad_n.id, nt_pad_p.id, μd=0.0, χ=0.6, v_tol=v_tol)

finalize!(mech_scen)

###############

q_tra_z_0 = 0.10
q_rev_y_0 = pi/4
set_configuration!(mech_scen, nt_tra_z.joint, [q_tra_z_0])
set_configuration!(mech_scen, nt_rev_y.joint, [q_rev_y_0])
set_state_spq!(mech_scen, nt_penci.joint, trans=SVector(penci_length * -0.8, 0.0, 0.0), rot=RotZ(-pi/2))

###############

x = get_state(mech_scen)

### Add meshcat visualizer
if !@isdefined vis
    vis = Visualizer()
    open(vis)
end

if true  ### visualization
	color_gray  = RGBA{Float32}(0.5, 0.5, 0.5, 1.0)
	color_red   = RGBA{Float32}(1.0, 0.0, 0.0, 1.0)
	color_green = RGBA{Float32}(0.0, 1.0, 0.0, 1.0)
	color_blue  = RGBA{Float32}(0.0, 0.0, 1.0, 0.5)
	mvis = MechanismVisualizer(my_mechanism, vis)

	### Add meshes to visualizer
	set_body_mesh_visual!(mvis, mech_scen, nt_pad_n.id,    color_gray)
	set_body_mesh_visual!(mvis, mech_scen, nt_pad_p.id,    color_gray)
	set_body_mesh_visual!(mvis, mech_scen, nt_rev_y_id.id, color_gray)
	set_body_mesh_visual!(mvis, mech_scen, nt_penci.id,    color_blue)

	set_mesh_visual!(     mvis, mech_scen, nt_plane.id,    color_green)
	set_configuration!(mech_scen, mvis)
end

### Run forward dynamics
rr = Radau_for_MechanismScenario(mech_scen)
rr.step.h_min = 1.0e-14
rr.step.h_max = 0.01
rr.step.h = rr.step.h_max

t_final   = 18.0
max_steps = 8200
println()

@time data_time, data_state = integrate_scenario_radau(rr, t_final=t_final, max_steps=max_steps)
println("t_sim: ", data_time[end])
println("simulated steps: ", length(data_time))
println("last time steps: ", (data_time[end] - data_time[end-10]) / 10)

### Move camera
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 7)
settransform!(vis["/Cameras/default"], Translation(0.0, 0.0, 0.30) ∘ LinearMap(RotZ(-π * 0.35)))

### Playback data
play_recorded_data(mvis, mech_scen, data_time, data_state, slowdown=1.0)
