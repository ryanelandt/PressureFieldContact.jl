@RigidBodyDynamics.indextype BristleID

abstract type FrictionModel end

struct Bristle <: FrictionModel
    BristleID::BristleID
    τ::Float64
    k̄::Float64
    Bristle(bristle_ID::BristleID; τ::Float64, k̄::Float64) = new(bristle_ID, τ, k̄)
end

struct Regularized <: FrictionModel
    v_tol⁻¹::Float64
    Regularized(v_tol) = new(1 / v_tol)
end

struct ContactInstructions
    id_1::MeshID
    id_2::MeshID
    μ::Float64
    χ::Float64
    FrictionModel::Union{Regularized,Bristle}
    function ContactInstructions(id_tri::MeshID, id_tet::MeshID, fric_model::Union{Regularized,Bristle};
            μ::Float64, χ::Float64)

        (0.0 <= μ <= 3.0) || error("mu in unexpected range.")
        return new(id_tri, id_tet, μ, χ, fric_model)
    end
end

struct TractionCache{N,T}
    n̂::FreeVector3D{SVector{3,T}}
    r_cart::NTuple{N,Point3D{SVector{3,T}}}
    v_cart::NTuple{N,FreeVector3D{SVector{3,T}}}
    dA::NTuple{N,T}
    p::NTuple{N,T}
end

@inline calc_p_dA(t::TractionCache, k::Int64) = t.p[k] * t.dA[k]

mutable struct spatialStiffness{T}
    K::Hermitian{T,Matrix{T}}
    σ_sqrt::Diagonal{T,MVector{6,T}}
    K⁻¹_sqrt::MMatrix{6,6,T,36}
    mul_pre::MMatrix{6,6,T,36}
    function spatialStiffness{T}() where {T}
        K = Hermitian(Matrix(T.(zeros(6,6) .+ NaN)))
        σ_sqrt = Diagonal(MVector{6,T}(NaN, NaN, NaN, NaN, NaN, NaN))
        K⁻¹_sqrt = MMatrix{6,6,T,36}(T.(zeros(6,6) .+ NaN))
        mul_pre = MMatrix{6,6,T,36}(T.(zeros(6,6)) .+ NaN)
        return new(K, σ_sqrt, K⁻¹_sqrt, mul_pre)
    end
end

mutable struct TypedElasticBodyBodyCache{N,T}
    spatialStiffness::spatialStiffness{T}
    quad::TriTetQuadRule{3,N}
    TractionCache::VectorCache{TractionCache{N,T}}
    mesh_1::MeshCache
    mesh_2::MeshCache
    x_rʷ_r²::Transform3D{T}
    x_r²_rʷ::Transform3D{T}
    x_r¹_r²::Transform3D{T}
	x_r²_r¹::Transform3D{T}
    twist_r²_r¹_r²::Twist{T}
    χ::Float64
    μ::Float64
    Ē::Float64
    wrench::Wrench{T}
    function TypedElasticBodyBodyCache{N,T}(quad::TriTetQuadRule{3,N}) where {N,T}
        trac_cache = VectorCache{TractionCache{N, T}}()
        K = spatialStiffness{T}()
        return new{N,T}(K, quad, trac_cache)
    end
end

struct TypedMechanismScenario{N,T}
    state::MechanismState{T}
    s::SegmentedVector{BristleID,T,Base.OneTo{BristleID},Array{T,1}}
    result::DynamicsResult{T}
    ṡ::SegmentedVector{BristleID,T,Base.OneTo{BristleID},Array{T,1}}
    f_generalized::Vector{T}
    bodyBodyCache::TypedElasticBodyBodyCache{N,T}
    GeometricJacobian::RigidBodyDynamics.CustomCollections.CacheIndexDict{BodyID,Base.OneTo{BodyID},Union{Nothing,GeometricJacobian{Array{T,2}}}}
    torque_third_law::Vector{T}
    τ_ext::Vector{T}
    rhs::Vector{T}
    nv::Int64

    function TypedMechanismScenario{N,T}(mechanism::Mechanism, quad::TriTetQuadRule{3,N}, v_path, body_ids,
            n_bristle_pairs::Int64) where {N,T}

        function makeJacobian(v_path, state::MechanismState{T}, body_ids::Base.OneTo{BodyID}) where {T}
            v_jac = RigidBodyDynamics.BodyCacheDict{Union{Nothing,GeometricJacobian{Array{T,2}}}}(body_ids)
            fill_with_nothing!(v_jac)
            for (body_id_k, path_k) = v_path
                if path_k != nothing  # body has no meshes attached to it
                    if body_id_k != BodyID(root_body(state.mechanism))
                        new_jac = geometric_jacobian(state, path_k)
                        fill!(new_jac.linear, NaN + zero(T))
                        fill!(new_jac.angular, NaN + zero(T))
                        v_jac[body_id_k] = new_jac
                    end
                end
            end
            return v_jac
        end
        function_Int64_six(k) = 6

        state = MechanismState{T}(mechanism)
        result = DynamicsResult{T}(mechanism)
        f_generalized = Vector{T}(undef, num_positions(mechanism))
        bodyBodyCache = TypedElasticBodyBodyCache{N,T}(quad)
        v_jac = makeJacobian(v_path, state, body_ids)
        n_dof_bristle = 6 * n_bristle_pairs
        s = SegmentedVector{BristleID}(zeros(T, n_dof_bristle), Base.OneTo(BristleID(n_bristle_pairs)), function_Int64_six)
        ṡ = SegmentedVector{BristleID}(zeros(T, n_dof_bristle), Base.OneTo(BristleID(n_bristle_pairs)), function_Int64_six)
        torque_third_law = Vector{T}(undef, num_positions(mechanism))
        τ_ext = zeros(T, num_positions(mechanism))
        rhs = zeros(T, num_positions(mechanism))
        nv = num_velocities(mechanism)
        return new{N,T}(state, s, result, ṡ, f_generalized, bodyBodyCache, v_jac, torque_third_law, τ_ext, rhs, nv)
    end

	function TypedMechanismScenario{N,T}(mechanism::Mechanism) where {N,T}
		return new{N,T}(MechanismState{T}(mechanism))
	end
end

mutable struct DiscreteControl
	dt::Float64
	t_last::Float64
	control::Function
	# DiscreteControl(control::Nothing) = nothing
	DiscreteControl(control::Function, dt=0.01) = new(dt, 0.0, control)
end

mutable struct MechanismScenario{NQ,T}
	float::TypedMechanismScenario{NQ,Float64}
	dual::TypedMechanismScenario{NQ,T}
    body_ids::Base.OneTo{BodyID}
    mesh_ids::Base.OneTo{MeshID}
    bristle_ids::Base.OneTo{BristleID}
	MeshCache::RigidBodyDynamics.CustomCollections.CacheIndexDict{MeshID,Base.OneTo{MeshID},MeshCache}
	ContactInstructions::Vector{ContactInstructions}
	TT_Cache::TT_Cache
	de::Function
	continuous_controller::Union{Nothing,Function}
	discrete_controller::Union{Nothing,DiscreteControl}
	# discrete_controller::Union{Nothing,Function}

    τ_ext::Vector{Float64}
    path::RigidBodyDynamics.CustomCollections.IndexDict{BodyID,Base.OneTo{BodyID},Union{Nothing,RigidBodyDynamics.Graphs.TreePath{RigidBody{Float64},Joint{Float64,JT} where JT<:JointType{Float64}}}}

	function MechanismScenario(; de::Function=calcXd!, n_quad_rule::Int64=2, N_chunk::Int64=6,
			continuous_controller::Union{Nothing,Function}=nothing,
			discrete_controller::Union{Nothing,DiscreteControl}=nothing,
			gravity::SVector{3,Float64}=SVector{3,Float64}(0.0, 0.0, -round(9.8054, digits=8, base=2) ))
		#
		mechanism = Mechanism(RigidBody{Float64}("world"); gravity=gravity)  # create empty mechanism
        (1 <= n_quad_rule <= 2) || error("only quadrature rules 1 (first order) and 2 (second? order) are currently implemented")
        quad = getTriQuadRule(n_quad_rule)
        NQ = length(quad.w)
        T = Dual{Nothing,Float64,N_chunk}
		cache_float = TypedMechanismScenario{NQ,Float64}(mechanism)
		cache_dual = TypedMechanismScenario{NQ,T}(mechanism)
		body_ids = Base.OneTo(BodyID(1))
		mesh_ids = Base.OneTo(MeshID(0))
		bristle_ids = Base.OneTo(BristleID(0))
		mesh_cache = MeshCacheDict{MeshCache}(mesh_ids)
		c_ins = Vector{ContactInstructions}()

		return new{NQ,T}(cache_float, cache_dual, body_ids, mesh_ids, bristle_ids, mesh_cache, c_ins, TT_Cache(), de,
			continuous_controller, discrete_controller)
    end
end

function finalize!(m::MechanismScenario{NQ,T}) where {NQ,T}
	function makePaths(mechanism::Mechanism, mesh_cache::MeshCacheDict{MeshCache}, body_ids::Base.OneTo{BodyID})
		the_type = Union{Nothing,RigidBodyDynamics.Graphs.TreePath{RigidBody{Float64},Joint{Float64,JT} where JT<:JointType{Float64}}}
		v_path = RigidBodyDynamics.BodyDict{the_type}(body_ids)
		fill_with_nothing!(v_path)
		for mesh_k = values(mesh_cache)
			body_id = mesh_k.BodyID
			if body_id != BodyID(root_body(mechanism))
				v_path[body_id] = path(mechanism, root_body(mechanism), bodies(mechanism)[body_id])
			end
		end
		return v_path
	end

	quad = getTriQuadRule(ifelse(NQ == 1, 1, 2))
	mechanism = get_mechanism(m)
	body_ids = Base.OneTo(last(bodies(mechanism)).id)
	m.body_ids = body_ids
	n_bristle_pairs = length(m.bristle_ids)
	m.τ_ext = zeros(Float64, num_positions(mechanism))
	m.path = makePaths(mechanism, m.MeshCache, body_ids)
	cache_path = m.path
	m.float = TypedMechanismScenario{NQ,Float64}(mechanism, quad, cache_path, body_ids, n_bristle_pairs)
	m.dual  = TypedMechanismScenario{NQ,T}(mechanism, quad, cache_path, body_ids, n_bristle_pairs)
	return nothing
end

num_partials(m::MechanismScenario{NQ,Dual{Nothing,Float64,N_partials}}) where {NQ,N_partials} = N_partials
function num_x(m::MechanismScenario{NQ,T}) where {NQ,T}
    tm = m.float
    mechanism = tm.state.mechanism
    return num_positions(mechanism) + num_velocities(mechanism) + length(tm.s)
 end
type_dual(m::MechanismScenario{NQ,T}) where {NQ,T} = T

function get_state(m::MechanismScenario)
    x = zeros(num_x(m))
    copyto!(x, m.float)
    return x
end

function set_state_spq!(m::MechanismScenario, joint::Joint; rot::Rotation=one(Quat{Float64}),
        trans::SVector{3,Float64}=zeros(SVector{3,Float64}), w::SVector{3,Float64}=zeros(SVector{3,Float64}),
        vel::SVector{3,Float64}=zeros(SVector{3,Float64}))

    rot = components(SPQuat(rot))
    state = m.float.state
    set_configuration!(state, joint, vcat(rot, trans))
    set_velocity!(state, joint, vcat(w, vel))
    return nothing
end

function addMesh!(m::MechanismScenario, mesh::MeshCache)
    mesh_ids_old = m.mesh_ids
    mesh_ids_new = Base.OneTo(MeshID(length(m.mesh_ids) + 1))
    mesh_cache = MeshCacheDict{MeshCache}(mesh_ids_new)
    for id = mesh_ids_old
        mesh_cache[id] = m.MeshCache[id]
    end
    mesh_cache[mesh_ids_new[end]] = mesh
    m.MeshCache = mesh_cache
    m.mesh_ids = mesh_ids_new
    return nothing
end

get_mechanism(m::MechanismScenario) = m.float.state.mechanism

function add_body_contact!(m::MechanismScenario, name::String, e_mesh::eMesh;
        i_prop::InertiaProperties,
        c_prop::Union{Nothing,ContactProperties}=nothing,
        body::Union{RigidBody{Float64},Nothing}=nothing,
        joint::JT=SPQuatFloating{Float64}(),
        dh::basic_dh=one(basic_dh{Float64})) where {JT<:JointType}

    nt = add_body!(m, name, e_mesh, i_prop=i_prop, body=body, joint=joint, dh=dh)
    nt_new_contact = add_contact!(m, name, e_mesh, c_prop=c_prop, body=nt.body)
    return NamedTuple{(:body, :joint, :id)}((nt.body, nt.joint, nt_new_contact.id))
end

# function make_eTree_obb(eM_box::eMesh{T1,T2}, c_prop::Union{Nothing,ContactProperties}) where {T1,T2}
#     xor(c_prop == nothing, T2 == Nothing) && error("Attempting to use nothing as the ContartProperties for a Tet mesh")
#     e_tree = eTree(eM_box, c_prop)
#     if T1 != Nothing
#         all_obb = [fit_tri_obb(eM_box, k) for k = 1:n_tri(eM_box)]
#     else
# 		all_obb = [fit_tet_obb(eM_box, k) for k = 1:n_tet(eM_box)]
#     end
# 	tree = obb_tree_from_aabb(get_tree(e_tree), all_obb)
#     return eTree(tree, c_prop)
# end

function add_contact!(m::MechanismScenario, name::String, e_mesh::eMesh;  c_prop::Union{Nothing,ContactProperties}=nothing,
        body::Union{RigidBody{Float64},Nothing}=nothing)

	function verify_eMesh_ContactProperties(eM::eMesh{T1,T2}, c_prop::Union{Nothing,ContactProperties}) where {T1,T2}
		(T1 == Tri) && (T2 == Tet) && error("eMesh has triangles and tets. Use as_tri_eMesh or as_tet_eMesh to convert eMesh.")
		(T1 == nothing) && (T2 == nothing) && error("eMesh has neither triangles and tets.")
		(T1 == Tri) && (c_prop != nothing) && error("Using ContactProperties for triangular eMesh")
		(T2 == Tet) && (c_prop == nothing) && error("Using nothing for tet eMesh")
	end

    verify_eMesh_ContactProperties(e_mesh, c_prop)
    body = return_body_never_nothing(get_mechanism(m), body)
    # e_tree = make_eTree_obb(e_mesh, c_prop)
	e_tree = eTree(e_mesh, c_prop)
    mesh = MeshCache(name, e_mesh, e_tree, body)
    addMesh!(m, mesh)
    return NamedTuple{(:id,)}((find_mesh_id(m, mesh),))
end

function add_body!(m::MechanismScenario, name::String, e_mesh::eMesh; i_prop::InertiaProperties,
        body::Union{RigidBody{Float64},Nothing}=nothing, joint::JT=SPQuatFloating{Float64}(),
        dh::basic_dh=one(basic_dh{Float64})) where {JT<:JointType}

    mesh_inertia_info = makeInertiaInfo(e_mesh, i_prop)
    return add_body_from_inertia!(get_mechanism(m), name, mesh_inertia_info, joint=joint, body=body, dh=dh)
end

return_body_never_nothing(mechanism::Mechanism, body::Nothing) = root_body(mechanism)
return_body_never_nothing(mechanism::Mechanism, body::RigidBody{Float64}) = body

function add_body_from_inertia!(mechanism::Mechanism, name::String, mesh_inertia_info::MeshInertiaInfo;
        joint::JT=SPQuatFloating{Float64}(), body::Union{RigidBody{Float64},Nothing}=nothing,
        dh::basic_dh{Float64}=one(basic_dh{Float64})) where {JT<:JointType}
    # TODO: check that a spherical floating joint isn't added

    body_parent = return_body_never_nothing(mechanism, body)
    body_child = newBodyFromInertia(name, mesh_inertia_info)
    j_parent_child, x_parent_child = outputJointTransform_ParentChild(body_parent, body_child, joint, dh)
    attach!(mechanism, body_parent, body_child, j_parent_child, joint_pose=x_parent_child)
    return NamedTuple{(:body, :joint)}((body_child, j_parent_child))
end

default_χ() = 0.5
default_μ() = 0.3

function add_friction_regularize!(m::MechanismScenario, mesh_id_1::MeshID, mesh_id_2::MeshID; μ::Float64=default_μ(),
		χ::Float64=default_χ(), v_tol::Float64=0.01)

    regularized = Regularized(v_tol)
    return add_friction!(m, mesh_id_1, mesh_id_2, regularized, μ=μ, χ=χ)
end

function add_friction_bristle!(m::MechanismScenario, mesh_id_1::MeshID, mesh_id_c::MeshID; τ::Float64=0.05, k̄=1.0e4,
		μ::Float64=default_μ(), χ::Float64=default_χ())

    isa(μ, Nothing) || (0 < μ) || error("μ cannot be 0 for bristle friction")
    bristle_id = BristleID(1 + length(m.bristle_ids))
    bf = Bristle(bristle_id, τ=τ, k̄=k̄)
    m.bristle_ids = Base.OneTo(bristle_id)
    return add_friction!(m, mesh_id_1, mesh_id_c, bf, μ=μ, χ=χ)
end

add_friction!(m::MechanismScenario, id_1::MeshID, id_2::MeshID, fric_model::FrictionModel; μ::Float64,
	χ::Float64) = add_friction!(m, id_1, id_2, m.MeshCache[id_1], m.MeshCache[id_2], fric_model, μ=μ, χ=χ)

function add_friction!(m::MechanismScenario, id_1::MeshID, id_2::MeshID, m_1::MeshCache{Nothing,Tet},
	m_2::MeshCache{Tri,Nothing}, fric_model::Union{Regularized,Bristle}; μ::Float64, χ::Float64)

	return add_friction!(m, id_2, id_1, fric_model; μ=μ, χ=χ)
end

function add_friction!(m::MechanismScenario, id_1::MeshID, id_2::MeshID, m_1::MeshCache{T1,T2},
	m_2::MeshCache{Nothing,Tet}, fric_model::Union{Regularized,Bristle}; μ::Float64, χ::Float64) where {T1,T2}

	push!(m.ContactInstructions, ContactInstructions(id_1, id_2, fric_model, μ=μ, χ=χ))
end
