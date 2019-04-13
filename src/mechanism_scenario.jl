
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
    frame_world::CartesianFrame3D
    quad::TriTetQuadRule{3,N}
    TractionCache::VectorCache{TractionCache{N,T}}
    mesh_1::MeshCache
    mesh_2::MeshCache
    x_rʷ_r¹::Transform3D{T}
    x_rʷ_r²::Transform3D{T}
    x_r²_rʷ::Transform3D{T}
    x_r¹_rʷ::Transform3D{T}
    x_r¹_r²::Transform3D{T}
	x_r²_r¹::Transform3D{T}
    twist_r²_r¹::Twist{T}
    twist_r²_r¹_r²::Twist{T}
    χ::Float64
    μ::Float64
    Ē::Float64
    wrench::Wrench{T}
    function TypedElasticBodyBodyCache{N,T}(frame_world::CartesianFrame3D, quad::TriTetQuadRule{3,N}) where {N,T}
        trac_cache = VectorCache{TractionCache{N, T}}()
        K = spatialStiffness{T}()
        return new{N,T}(K, frame_world, quad, trac_cache)
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
    frame_world::CartesianFrame3D  # TODO: does this need to be here?
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
        frame_world = root_frame(mechanism)
        bodyBodyCache = TypedElasticBodyBodyCache{N,T}(frame_world, quad)
        v_jac = makeJacobian(v_path, state, body_ids)
        n_dof_bristle = 6 * n_bristle_pairs
        s = SegmentedVector{BristleID}(zeros(T, n_dof_bristle), Base.OneTo(BristleID(n_bristle_pairs)), function_Int64_six)
        ṡ = SegmentedVector{BristleID}(zeros(T, n_dof_bristle), Base.OneTo(BristleID(n_bristle_pairs)), function_Int64_six)
        torque_third_law = Vector{T}(undef, num_positions(mechanism))
        τ_ext = zeros(T, num_positions(mechanism))
        rhs = zeros(T, num_positions(mechanism))
        nv = num_velocities(mechanism)
        return new{N,T}(state, s, result, ṡ, f_generalized, bodyBodyCache, v_jac, root_frame(mechanism), torque_third_law, τ_ext, rhs, nv)
    end
end

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

struct MechanismScenario{NX,NQ,T}
    body_ids::Base.OneTo{BodyID}
    mesh_ids::Base.OneTo{MeshID}
    bristle_ids::Base.OneTo{BristleID}
    frame_world::CartesianFrame3D
    TT_Cache::TT_Cache
    τ_ext::Vector{Float64}
    float::TypedMechanismScenario{NQ,Float64}
    dual::TypedMechanismScenario{NQ,T}
    path::RigidBodyDynamics.CustomCollections.IndexDict{BodyID,Base.OneTo{BodyID},Union{Nothing,RigidBodyDynamics.Graphs.TreePath{RigidBody{Float64},Joint{Float64,JT} where JT<:JointType{Float64}}}}
    MeshCache::RigidBodyDynamics.CustomCollections.CacheIndexDict{MeshID,Base.OneTo{MeshID},MeshCache}
    ContactInstructions::Vector{ContactInstructions}
    de::Function
    continuous_controller::Union{Nothing,Function}
    discrete_controller::Union{Nothing,Function}
    function MechanismScenario(ts::TempContactStruct, de::Function; continuous_controller::Union{Nothing,Function}=nothing,
            discrete_controller::Union{Nothing,Function}=nothing, n_quad_rule::Int64=2, N_chunk::Int64=6)

        (1 <= n_quad_rule <= 2) || error("only quadrature rules 1 (first order) and 2 (second? order) are currently implemented")
        quad = getTriQuadRule(n_quad_rule)
        NQ = length(quad.w)

        mechanism = ts.mechanism
        body_ids = Base.OneTo(last(bodies(mechanism)).id)
        bristle_ids = ts.bristle_ids
        n_bristle_pairs = length(bristle_ids)

        n_q = num_positions(mechanism)
        n_v = num_velocities(mechanism)
        (n_q == n_v) || error("n_q ($n_q) and n_v ($n_v) are different. Replace QuaternionFloating joints with SPQuatFloating joints.")

        NX = n_q + n_v + 6 * n_bristle_pairs

        T = Dual{Nothing,Float64,N_chunk}
        frame_world = root_frame(mechanism)
        τ_ext = zeros(Float64, num_positions(mechanism))
        mesh_cache = ts.MeshCache
        cache_path = makePaths(mechanism, mesh_cache, body_ids)
        cache_float = TypedMechanismScenario{NQ,Float64}(mechanism, quad, cache_path, body_ids, n_bristle_pairs)
        cache_dual = TypedMechanismScenario{NQ,T}(mechanism, quad, cache_path, body_ids, n_bristle_pairs)
        vec_instructions = ts.ContactInstructions
        return new{NX,NQ,T}(body_ids, ts.mesh_ids, bristle_ids, frame_world, TT_Cache(), τ_ext, cache_float, cache_dual,
            cache_path, mesh_cache, vec_instructions, de, continuous_controller, discrete_controller)  # MVector{1,Float64}(0.0))
    end
end

num_partials(m::MechanismScenario{NX,NQ,Dual{Nothing,Float64,N_partials}}) where {NX,NQ,N_partials} = N_partials
num_x(m::MechanismScenario{NX,NQ,T}) where {NX,NQ,T} = NX
type_dual(m::MechanismScenario{NX,NQ,T}) where {NX,NQ,T} = T

function get_state(mech_scen::MechanismScenario)
    x = zeros(num_x(mech_scen))
    copyto!(x, mech_scen.float)
    return x
end

function set_state_spq!(mech_scen::MechanismScenario, joint::Joint; rot::Rotation=one(Quat{Float64}),
        trans::SVector{3,Float64}=zeros(SVector{3,Float64}), w::SVector{3,Float64}=zeros(SVector{3,Float64}),
        vel::SVector{3,Float64}=zeros(SVector{3,Float64}))

    rot = components(SPQuat(rot))
    state = mech_scen.float.state
    set_configuration!(state, joint, vcat(rot, trans))
    set_velocity!(state, joint, vcat(w, vel))
    return nothing
end
