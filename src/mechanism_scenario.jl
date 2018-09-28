
@RigidBodyDynamics.indextype BristleID

struct BristleFriction
    BristleID::BristleID
    tau::Float64
    K_r::Float64
    K_theta::Float64
    function BristleFriction(bristle_ID::BristleID, tau::Float64, K_r::Float64, K_theta::Float64)
        return new(bristle_ID, tau, K_r, K_theta)
    end
end

mutable struct ContactInstructions
    id_tri::MeshID
    id_tet::MeshID
    frac_epsilon::Float64
    frac_linear_weight::Float64
    mu_pair::Float64
    BristleFriction::Union{Nothing,BristleFriction}
    function ContactInstructions(id_tri::MeshID, id_tet::MeshID, frac_epsilon::Float64, frac_linear_weight::Float64, mu_pair::Float64)
        return new(id_tri, id_tet, frac_epsilon, frac_linear_weight, mu_pair, nothing)
    end
    function ContactInstructions(id_tri::MeshID, id_tet::MeshID, frac_epsilon::Float64, frac_linear_weight::Float64, mu_pair::Float64, bristle_friction::BristleFriction)
        return new(id_tri, id_tet, frac_epsilon, frac_linear_weight, mu_pair, bristle_friction)
    end
end

mutable struct TractionCache{N,T}
    traction_normal::FreeVector3D{SVector{3,T}}
    r_cart::MVector{N,Point3D{SVector{3,T}}}
    v_cart_t::MVector{N,FreeVector3D{SVector{3,T}}}
    p_dA::MVector{N,T}
    function TractionCache{N,T}(frame::CartesianFrame3D) where {N,T}
        traction_normal = FreeVector3D(frame, SVector{3, T}(NaN .+ zeros(3)))
        r_cart = MVector{N,Point3D{SVector{3,Float64}}}([Point3D(frame, SVector{3,Float64}(NaN, NaN, NaN)) for k = 1:3])
        v_cart_t = MVector{N,FreeVector3D{SVector{3,Float64}}}([FreeVector3D(frame, SVector{3,Float64}(NaN, NaN, NaN)) for k = 1:3])
        p_dA = MVector{N, T}(NaN .+ zeros(N))
        return new(traction_normal, r_cart, v_cart_t, p_dA)
    end
end

mutable struct TypedQuadTriCache{T}
    clip_poly_4D_1::ClippedPolygon{4,T}
    clip_poly_4D_2::ClippedPolygon{4,T}
    clip_poly_3D::ClippedPolygon{3,T}
    area_quad_k::T
    A_zeta_phi::MatrixTransform{4,3,T,12}
    function TypedQuadTriCache{T}(frame_world::CartesianFrame3D) where {T}
        clip_poly_4D_1 = ClippedPolygon{4,T}(frame_tet_cs)
        clip_poly_4D_2 = ClippedPolygon{4,T}(frame_tet_cs)
        clip_poly_3D = ClippedPolygon{3,T}(frame_world)
        return new(clip_poly_4D_1, clip_poly_4D_2, clip_poly_3D)
    end
end

mutable struct TypedTriTetCache{T}
    quadTriCache::TypedQuadTriCache{T}
    strain::SVector{4,Float64}
    A_w_zeta_top::MatrixTransform{3,4,T,12}
    traction_normal::FreeVector3D{SVector{3,T}}
    centroid_zeta::Point4D{SVector{4,T}}
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
    frac_epsilon::Float64
    hc_velocity_damping::Float64
    mu::Float64
    E_effective::Float64
    inv_thickness::Float64
    function TypedElasticBodyBodyCache{N,T}(frame_world::CartesianFrame3D, quad::TriTetQuadRule{3,N}) where {N,T}
        triTetCache = TypedTriTetCache{T}(frame_world)
        trac_cache = VectorCache(TractionCache{N,T}(frame_world))
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
        s = SegmentedVector{BristleID}(Vector{T}(undef, n_dof_bristle), Base.OneTo(BristleID(n_bristle_pairs)), function_Int64_six)
        ṡ = SegmentedVector{BristleID}(Vector{T}(undef, n_dof_bristle), Base.OneTo(BristleID(n_bristle_pairs)), function_Int64_six)
        return new{N,T}(state, s, result, ṡ, f_generalized, bodyBodyCache, v_jac)
    end
end

struct MechanismScenario{N,T}
    body_ids::Base.OneTo{BodyID}
    mesh_ids::Base.OneTo{MeshID}
    bristle_ids::Base.OneTo{BristleID}
    frame_world::CartesianFrame3D
    TT_Cache::TT_Cache
    tau_ext::Vector{Float64}
    float::TypedMechanismScenario{N,Float64}
    dual::TypedMechanismScenario{N,T}
    path::RigidBodyDynamics.CustomCollections.IndexDict{BodyID,Base.OneTo{BodyID},Union{Nothing,RigidBodyDynamics.Graphs.TreePath{RigidBody{Float64},Joint{Float64,JT} where JT<:JointType{Float64}}}}
    MeshCache::RigidBodyDynamics.CustomCollections.CacheIndexDict{MeshID,Base.OneTo{MeshID},MeshCache}
    ContactInstructions::Vector{ContactInstructions}
    de::Function
    function MechanismScenario(de::Function, mechanism::Mechanism, vec_MeshCache::Vector{MeshCache}, T, n_bristle_pairs::Int64)
        function makeMeshCacheDict(mechanism::Mechanism, vec_MeshCache::Vector{MeshCache}, mesh_ids::Base.OneTo{MeshID})
            cache_mesh = MeshCacheDict{MeshCache}(mesh_ids)
            for id = mesh_ids
                cache_mesh[id] = vec_MeshCache[id]  # should index vec_MeshCache fine
            end
            return cache_mesh
        end
        function makePaths(mechanism::Mechanism, vec_MeshCache::Vector{MeshCache}, body_ids::Base.OneTo{BodyID})
            the_type = Union{Nothing,RigidBodyDynamics.Graphs.TreePath{RigidBody{Float64},Joint{Float64,JT} where JT<:JointType{Float64}}}
            v_path = RigidBodyDynamics.BodyDict{the_type}(body_ids)
            fill_with_nothing!(v_path)
            for mesh_k = vec_MeshCache
                body_id = mesh_k.BodyID
                if body_id != BodyID(root_body(mechanism))
                    v_path[body_id] = path(mechanism, root_body(mechanism), bodies(mechanism)[body_id])
                end
            end
            return v_path
        end

        quad = getTriQuadRule(2)  # TODO: move somewhere else
        N = length(quad.w)
        tau_ext = zeros(Float64, num_positions(mechanism))
        body_ids = Base.OneTo(BodyID(num_bodies(mechanism)))
        mesh_ids = Base.OneTo(MeshID(length(vec_MeshCache)))
        cache_mesh = makeMeshCacheDict(mechanism, vec_MeshCache, mesh_ids)
        cache_path = makePaths(mechanism, vec_MeshCache, body_ids)
        vec_instructions = Vector{ContactInstructions}()
        cache_float = TypedMechanismScenario{N,Float64}(mechanism, quad, cache_path, body_ids, n_bristle_pairs)
        cache_dual = TypedMechanismScenario{N,T}(mechanism, quad, cache_path, body_ids, n_bristle_pairs)
        frame_world = root_frame(mechanism)
        bristle_ids = Base.OneTo(n_bristle_pairs)
        return new{N,T}(body_ids, mesh_ids, bristle_ids, frame_world, TT_Cache(), tau_ext, cache_float, cache_dual, cache_path, cache_mesh, vec_instructions, de)
    end
end
