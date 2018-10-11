
struct TractionCache{N,T}
    traction_normal::FreeVector3D{SVector{3,T}}
    r_cart::NTuple{N,Point3D{SVector{3,T}}}
    v_cart_t::NTuple{N,FreeVector3D{SVector{3,T}}}
    p_dA::NTuple{N,T}
    function TractionCache(traction_normal::FreeVector3D{SVector{3,T}}, r_cart::NTuple{N,Point3D{SVector{3,T}}},
                                v_cart_t::NTuple{N,FreeVector3D{SVector{3,T}}}, p_dA::NTuple{N,T}) where {N,T}

        return new{N,T}(traction_normal, r_cart, v_cart_t, p_dA)
    end
    function TractionCache(N, T)  # TODO: delete this when vector cache only works for immutables
        frame = frame_tet_cs

        traction_normal = FreeVector3D(frame, SVector{3, T}(NaN .+ zeros(3)))
        r_cart = NTuple{N,Point3D{SVector{3,T}}}([Point3D(frame, SVector{3,T}(NaN, NaN, NaN)) for k = 1:N])
        v_cart_t = NTuple{N,FreeVector3D{SVector{3,T}}}([FreeVector3D(frame, SVector{3,T}(NaN, NaN, NaN)) for k = 1:N])
        p_dA = NTuple{N,T}(NaN .+ zeros(N))
        return new{N,T}(traction_normal, r_cart, v_cart_t, p_dA)
    end
end

mutable struct TypedQuadTriCache{T}
    clip_poly_4D_1::ClippedPolygon{4,T}
    clip_poly_4D_2::ClippedPolygon{4,T}
    clip_poly_3D::ClippedPolygon{3,T}
    area_quad_k::T
    function TypedQuadTriCache{T}(frame_world::CartesianFrame3D) where {T}
        clip_poly_4D_1 = ClippedPolygon{4,T}(frame_tet_cs)
        clip_poly_4D_2 = ClippedPolygon{4,T}(frame_tet_cs)
        clip_poly_3D = ClippedPolygon{3,T}(frame_world)
        return new(clip_poly_4D_1, clip_poly_4D_2, clip_poly_3D)
    end
end

mutable struct TypedTriTetCache{T}
    quadTriCache::TypedQuadTriCache{T}
    ϵ::SVector{4,Float64}
    traction_normal::FreeVector3D{SVector{3,T}}
    centroid_ζ::Point4D{SVector{4,T}}
    centroid_w::Point3D{SVector{3,T}}
    TypedTriTetCache{T}(frame_world::CartesianFrame3D) where {T} = new(TypedQuadTriCache{T}(frame_world))
end

mutable struct TypedElasticBodyBodyCache{N,T}
    quad::TriTetQuadRule{3,N}
    triTetCache::TypedTriTetCache{T}
    TractionCache::VectorCache{TractionCache{N,T}}
    mesh_tri::MeshCache
    mesh_tet::MeshCache
    x_root_tri::Transform3D{T}
    x_root_tet::Transform3D{T}
    x_tet_root::Transform3D{T}
    twist_tri_tet::Twist{T}
    frac_ϵ::Float64
    χ::Float64
    μ::Float64
    Ē::Float64
    d⁻¹::Float64
    function TypedElasticBodyBodyCache{N,T}(frame_world::CartesianFrame3D, quad::TriTetQuadRule{3,N}) where {N,T}
        triTetCache = TypedTriTetCache{T}(frame_world)
        tc = TractionCache(N, T)
        trac_cache = VectorCache(tc)
        return new{N,T}(quad, triTetCache, trac_cache)
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
    function TypedMechanismScenario{N,T}(mechanism::Mechanism, quad::TriTetQuadRule{3,N}, v_path, body_ids, n_bristle_pairs::Int64) where {N,T}
        function makeJacobian(v_path, state::MechanismState{T}, body_ids::Base.OneTo{BodyID}) where {T}
            v_jac = RigidBodyDynamics.BodyCacheDict{Union{Nothing,GeometricJacobian{Array{T,2}}}}(body_ids)
            fill_with_nothing!(v_jac)
            for (body_id_k, path_k) = v_path
                if body_id_k != BodyID(root_body(state.mechanism))
                    new_jac = geometric_jacobian(state, path_k)
                    new_jac.linear .= NaN
                    new_jac.angular .= NaN
                    v_jac[body_id_k] = new_jac
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
        return new{N,T}(state, s, result, ṡ, f_generalized, bodyBodyCache, v_jac)
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

struct MechanismScenario{N,T}
    body_ids::Base.OneTo{BodyID}
    mesh_ids::Base.OneTo{MeshID}
    bristle_ids::Base.OneTo{BristleID}
    frame_world::CartesianFrame3D
    TT_Cache::TT_Cache
    τ_ext::Vector{Float64}
    float::TypedMechanismScenario{N,Float64}
    dual::TypedMechanismScenario{N,T}
    path::RigidBodyDynamics.CustomCollections.IndexDict{BodyID,Base.OneTo{BodyID},Union{Nothing,RigidBodyDynamics.Graphs.TreePath{RigidBody{Float64},Joint{Float64,JT} where JT<:JointType{Float64}}}}
    MeshCache::RigidBodyDynamics.CustomCollections.CacheIndexDict{MeshID,Base.OneTo{MeshID},MeshCache}
    ContactInstructions::Vector{ContactInstructions}
    de::Function
    function MechanismScenario(ts::TempContactStruct, de::Function; n_quad_rule::Int64=2)
        (1 <= n_quad_rule <= 2) || error("only quadrature rules 1 (first order) and 2 (second? order) are currently implemented")

        quad = getTriQuadRule(n_quad_rule)
        N = length(quad.w)
        mechanism = ts.mechanism
        body_ids = Base.OneTo(last(bodies(mechanism)).id)
        mesh_ids = ts.mesh_ids
        bristle_ids = ts.bristle_ids
        n_bristle_pairs = length(bristle_ids)

        n_positions = num_positions(mechanism)
        n_velocities = num_velocities(mechanism)
        (n_positions == n_velocities) || error("n_positions ($n_positions) and n_velocities ($n_velocities) are different. Replace QuaternionFloating joints with SPQuatFloating joints.")

        n_dof = n_positions + n_velocities + 6 * n_bristle_pairs
        T = Dual{Nothing,Float64,n_dof}
        frame_world = root_frame(mechanism)
        τ_ext = zeros(Float64, num_positions(mechanism))
        mesh_cache = ts.MeshCache
        cache_path = makePaths(mechanism, mesh_cache, body_ids)
        cache_float = TypedMechanismScenario{N,Float64}(mechanism, quad, cache_path, body_ids, n_bristle_pairs)
        cache_dual = TypedMechanismScenario{N,T}(mechanism, quad, cache_path, body_ids, n_bristle_pairs)
        vec_instructions = ts.ContactInstructions
        return new{N,T}(body_ids, mesh_ids, bristle_ids, frame_world, TT_Cache(), τ_ext, cache_float, cache_dual, cache_path, mesh_cache, vec_instructions, de)
    end
end

num_partials(m::MechanismScenario{N, Dual{Nothing,Float64,N_partials}}) where {N, N_partials} = N_partials

type_dual(m::MechanismScenario{N,T}) where {N,T} = T

function get_state(mech_scen::MechanismScenario)
    x = zeros(num_partials(mech_scen))
    copyto!(x, mech_scen.float)
    return x
end
