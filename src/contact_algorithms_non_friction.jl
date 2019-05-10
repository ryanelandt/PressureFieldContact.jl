
function calcXd!(xx::AbstractVector{T}, x::AbstractVector{T}, m::MechanismScenario{NQ,T}, t::Float64=0.0) where {NQ,T}
    return calcXd!(xx, x, m, m.dual, t)
end

function calcXd!(xx::AbstractVector{Float64}, x::AbstractVector{Float64}, m::MechanismScenario{NQ,T}, t::Float64=0.0) where {NQ,T}
    return calcXd!(xx, x, m, m.float, t)
end

"""
    calcXd!(xx, x, m, tm, t)

Conventions:
n̂ refers to the contact surface normal that points into body B
v_cart refers to + velocity of B - the velocity of A
the wrench is the wrench applied TO body A
"""
function calcXd!(xx::AbstractVector{T1}, x::AbstractVector{T1}, m::MechanismScenario{NQ,T2},
        tm::TypedMechanismScenario{NQ,T1}, t::Float64=0.0) where {NQ,T1,T2}

    state = tm.state
    copyto!(tm, x)
    H = tm.result.massmatrix
    mass_matrix!(H, state)
    dynamics_bias!(tm.result, state)
    configuration_derivative!(tm.result.q̇, state)
    forceAllElasticIntersections!(m, tm)

    (m.continuous_controller == nothing) || m.continuous_controller(tm, t)

	sum_all_forces!(m, tm)

	chol_fact = LinearAlgebra.cholesky!(H)
	ldiv!(tm.result.v̇.parent, chol_fact, tm.rhs)

    copyto!(xx, tm, tm.result)
    return nothing
end

function sum_all_forces!(m::MechanismScenario{NQ,T2}, tm::TypedMechanismScenario{NQ,Float64}) where {NQ,T2}
	BLAS.blascopy!(tm.nv, tm.f_generalized, 1, tm.rhs, 1)
	BLAS.axpy!(-1.0, tm.result.dynamicsbias.parent, tm.rhs)
	BLAS.axpy!(+1.0, tm.τ_ext, tm.rhs)
	BLAS.axpy!(+1.0, m.τ_ext, tm.rhs)
end

function sum_all_forces!(m::MechanismScenario{NQ,T2}, tm::TypedMechanismScenario{NQ,T1}) where {NQ,T1,T2}
	tm.rhs  .= tm.f_generalized
	tm.rhs .-= tm.result.dynamicsbias.parent
	tm.rhs .+= tm.τ_ext
	tm.rhs .+= m.τ_ext
end

function calcXd(x::AbstractVector{T1}, m::MechanismScenario{NQ,T2}, t::Float64=0.0) where {T1,NQ,T2}
    xx = deepcopy(x)
    calcXd!(xx, x, m)
    return xx
end

function forceAllElasticIntersections!(m::MechanismScenario{NQ,T1}, tm::TypedMechanismScenario{NQ,T2}) where {NQ,T1,T2}
    refreshJacobians!(m, tm)
    tm.f_generalized .= zero(T2)
    for k = 1:length(m.ContactInstructions)
        c_ins = m.ContactInstructions[k]
        force_single_elastic_intersection!(m, tm, c_ins)
    end
    return nothing
end

function force_single_elastic_intersection!(m::MechanismScenario{NQ,T1}, tm::TypedMechanismScenario{NQ,T2},
    c_ins::ContactInstructions) where {NQ,T1,T2}

    calcTriTetIntersections!(m, c_ins)
    if (0 != length(m.TT_Cache))  # yes intersections
        refreshBodyBodyCache!(m, tm, c_ins)
        integrate_over!(tm.bodyBodyCache, m.TT_Cache)
        if !isempty(tm.bodyBodyCache.TractionCache)
            wrench = yes_contact!(c_ins.FrictionModel, tm, c_ins)
            addGeneralizedForcesThirdLaw!(wrench, tm, c_ins)
            return wrench
        end
    end
    no_contact!(c_ins.FrictionModel, tm, c_ins)
end

function refreshJacobians!(m::MechanismScenario{NQ,T1}, tm::TypedMechanismScenario{NQ,T2}) where {NQ,T1,T2}
    for k = m.body_ids
        path_k = m.path[k]
        (path_k != nothing) && geometric_jacobian!(tm.GeometricJacobian[k], tm.state, path_k)
    end
    return nothing
end

function calcTriTetIntersections!(m::MechanismScenario, c_ins::ContactInstructions) # where {N,T}
    b = m.float.bodyBodyCache  # this can be float because intersection is assumed to not depend on partials
    refreshBodyBodyTransform!(m, m.float, c_ins)  # TODO: isn't b.x_r¹_rʷ already calculated?
	x_r¹_r² = b.x_r¹_r²
    update_TT_Cache!(m.TT_Cache, translation(x_r¹_r²), rotation(x_r¹_r²))
    tree_tree_intersect(m.TT_Cache, get_tree(b.mesh_1), get_tree(b.mesh_2))
    return nothing
end

function refreshBodyBodyTransform!(m::MechanismScenario, tm::TypedMechanismScenario{N,T},
        c_ins::ContactInstructions) where {N,T}

    b = tm.bodyBodyCache
    b.mesh_1 = m.MeshCache[c_ins.id_1]
    b.mesh_2 = m.MeshCache[c_ins.id_2]
      x_rʷ_r¹ = transform_to_root(tm.state, b.mesh_1.BodyID)  # TODO: add safe=false
    b.x_rʷ_r² = transform_to_root(tm.state, b.mesh_2.BodyID)  # TODO: add safe=false  # used for wrench calculation
    b.x_r²_rʷ = inv(b.x_rʷ_r²)
	b.x_r²_r¹ = b.x_r²_rʷ * x_rʷ_r¹
	b.x_r¹_r² = inv(b.x_r²_r¹)
    return nothing
end

function refreshBodyBodyCache!(m::MechanismScenario, tm::TypedMechanismScenario{N,T},
        c_ins::ContactInstructions) where {N,T}

    b = tm.bodyBodyCache
    empty!(b.TractionCache)
    refreshBodyBodyTransform!(m, tm, c_ins)

    # Rates
    twist_w_r¹ = twist_wrt_world(tm.state, b.mesh_1.BodyID)
    twist_w_r² = twist_wrt_world(tm.state, b.mesh_2.BodyID)
	  twist_r²_r¹ = -twist_w_r¹ + twist_w_r²  # velocity of tet wrt tri expressed in world
	b.twist_r²_r¹_r² = transform(twist_r²_r¹, b.x_r²_rʷ)

	b.μs = c_ins.μs
	b.μd = c_ins.μd
    b.χ = c_ins.χ
    b.Ē = get_c_prop(b.mesh_2).Ē
    return nothing
end

function integrate_over!(b::TypedElasticBodyBodyCache{N,T}, ttCache::TT_Cache) where {N,T}
	mesh_1 = b.mesh_1
    mesh_2 = b.mesh_2
	for k = 1:length(ttCache.vc)
		i_1, i_2 = ttCache.vc[k]
	    integrate_over!(i_1, i_2, mesh_1, mesh_2, b)
	end
end

function triangle_vertices(i_tri::Int64, m::MeshCache)
    ind_vert = get_ind_tri(m)[i_tri]
    return get_point(m)[ind_vert]
end

function tetrahedron_vertices_ϵ(i_tet::Int64, m::MeshCache)
    ind_vert = get_ind_tet(m)[i_tet]
    ϵ = get_ϵ(m)[ind_vert]
    ϵ = SMatrix{1,4,Float64,4}(ϵ)
    cart_vert = get_point(m)[ind_vert]
    return cart_vert, ϵ
end

function calc_ζ_transforms(p_tet)
	x_r¹_ζ = asMatOnePad(p_tet)
    x_ζ_r¹ = inv(x_r¹_ζ)
    return x_r¹_ζ, x_ζ_r¹
end

find_plane_tet(E::Float64, ϵ::SMatrix{1,4,Float64,4}, X_r_w) = (E * ϵ) * X_r_w

function integrate_over!(i_1::Int64, i_2::Int64, mesh_1::MeshCache{Nothing,Tet}, mesh_2::MeshCache{Nothing,Tet},
        b::TypedElasticBodyBodyCache{N,T}) where {N,T}

    vert_1, ϵ¹ = tetrahedron_vertices_ϵ(i_1, mesh_1)
    vert_2, ϵ² = tetrahedron_vertices_ϵ(i_2, mesh_2)
    x_r¹_ζ¹, x_ζ¹_r¹ = calc_ζ_transforms(vert_1)
    x_r²_ζ², x_ζ²_r² = calc_ζ_transforms(vert_2)

	ϵ_plane_1_r² = find_plane_tet(get_Ē(mesh_1), ϵ¹, x_ζ¹_r¹ * b.x_r¹_r².mat)  # TODO: is there a faster way to do this multiplication?
	ϵ_r² = ϵ² * x_ζ²_r²
    ϵ_plane_2_r² = find_plane_tet(get_Ē(mesh_2), ϵ², x_ζ²_r²)
	ϵ_plane_r² = ϵ_plane_2_r² - ϵ_plane_1_r²  # normalize penetration extent is positive so this describes the plane
		# of the contact surface pointing towards mesh_2

	x_r²_ζ¹ = b.x_r²_r¹.mat * x_r¹_ζ¹

    poly_r² = clip_plane_tet(ϵ_plane_r², x_r²_ζ¹)
    if 3 <= length(poly_r²)
        poly_ζ² = one_pad_then_mul(x_ζ²_r², poly_r²)
        poly_ζ² = zero_small_coordinates(poly_ζ²)  # This needs to be done to avoid a degenerate situation where the
            # plane lies exactly on the intersection of the faces of two tets. This situation happens EVERY time two
            # tet faces that lie on the surface intersect.
        poly_ζ² = clip_in_tet_coordinates(poly_ζ²)
        if 3 <= length(poly_ζ²)
			n̂² = normalize(unPad(ϵ_plane_r²))
			integrate_over_polygon_patch!(b, n̂², poly_ζ², x_r²_ζ², ϵ_r²)
        end
    end
end

function integrate_over!(i_1::Int64, i_2::Int64, mesh_1::MeshCache{Tri,Nothing}, mesh_2::MeshCache{Nothing,Tet},
        b::TypedElasticBodyBodyCache{N,T}) where {N,T}

    vert_1 = triangle_vertices(i_1, mesh_1)
    vert_2, ϵ² = tetrahedron_vertices_ϵ(i_2, mesh_2)
    x_r²_ζ², x_ζ²_r² = calc_ζ_transforms(vert_2)
	ϵ_r² = ϵ² * x_ζ²_r²

	x_ζ²_r¹ = x_ζ²_r² * b.x_r²_r¹.mat
	v1 = x_ζ²_r¹ * onePad(vert_1[1])
	v2 = x_ζ²_r¹ * onePad(vert_1[2])
	v3 = x_ζ²_r¹ * onePad(vert_1[3])
	poly_ζ² = clip_in_tet_coordinates(v1, v2, v3)

    if 3 <= length(poly_ζ²)
		n̂_r¹ = FreeVector3D(mesh_1.FrameID, triangleNormal(vert_1))
        n̂² = transform(n̂_r¹, b.x_r²_r¹)
		integrate_over_polygon_patch!(b, n̂².v, poly_ζ², x_r²_ζ², ϵ_r²)
    end
end

function integrate_over_polygon_patch!(b::TypedElasticBodyBodyCache{N,T}, n̂²::SVector{3,T},
		poly_ζ²::poly_eight{4,T},
        x_r²_ζ²::SMatrix{4,4,Float64,16},
		ϵ_r²::SMatrix{1,4,Float64,4}) where {N,T}

	poly_r² = mul_then_un_pad(x_r²_ζ², poly_ζ²)
	centroid_r² = centroid(poly_r², n̂²)[2]
    N_vertices = length(poly_ζ²)
	vert_r²_2 = poly_r²[N_vertices]
    for k = 1:N_vertices
		vert_r²_1 = vert_r²_2
		vert_r²_2 = poly_r²[k]
		area_quad_k = triangle_area((vert_r²_1, vert_r²_2, centroid_r²), n̂²)
		A_r²_ϕ = hcat(vert_r²_1, vert_r²_2, centroid_r²)
		(0.0 < area_quad_k) && fillTractionCacheForTriangle!(b, area_quad_k, n̂², A_r²_ϕ, ϵ_r²)  # no point in adding intersection if area is zero
    end
end

function fillTractionCacheForTriangle!(b::TypedElasticBodyBodyCache{N,T}, area_quad_k::T,
        n̂²::SVector{3,T},
		A_r²_ϕ::SMatrix{3,3,T,9},
		ϵ_r²::SMatrix{1,4,Float64,4}) where {N,T}

	for k = 1:N
		r_cart_1, dA_1, p_1 = fillTractionCacheInnerLoop!(k, b, area_quad_k, A_r²_ϕ, ϵ_r²)
		(0.0 < p_1) && addCacheItem!(b.TractionCache, TractionCache(n̂², r_cart_1, dA_1, p_1))
	end
end

@inline spatial_vel_formula(v::Twist{T}, b::SVector{3,T}) where {T} = linear(v) + cross(angular(v), b)

# TODO: write a function that performs "dot(ϵ_r², r²_pad.v)" without padding
function fillTractionCacheInnerLoop!(k::Int64, b::TypedElasticBodyBodyCache{N,T}, area_quad_k::T,
		A_r²_ϕ::SMatrix{3,3,T,9}, ϵ_r²::SMatrix{1,4,Float64,4}) where {N,T}

	rϕ = b.quad.zeta[k]
	r² = A_r²_ϕ * rϕ
	ϵ_quad = a_dot_one_pad_b(ϵ_r², r²)
	ṙ² = spatial_vel_formula(b.twist_r²_r¹_r², r²)
	ϵϵ = - a_dot_one_pad_b(ϵ_r², ṙ²)  # TODO: Why is there a negative sign?
	damp_term = max(zero(T), 1.0 + b.χ * ϵϵ)
    p_hc = ϵ_quad * b.Ē * damp_term
	dA = b.quad.w[k] * area_quad_k
	return r², dA, p_hc
end

function addGeneralizedForcesThirdLaw!(wrench::Wrench{T}, tm::TypedMechanismScenario{N,T}, cInfo::ContactInstructions) where {N,T}
    torque_third_law = tm.torque_third_law
    f_generalized = tm.f_generalized
	wrench = transform(wrench, tm.bodyBodyCache.x_rʷ_r²)
    addGeneralizedForcesExternal!(f_generalized,  wrench, tm.GeometricJacobian[tm.bodyBodyCache.mesh_1.BodyID], torque_third_law)
    addGeneralizedForcesExternal!(f_generalized, -wrench, tm.GeometricJacobian[tm.bodyBodyCache.mesh_2.BodyID], torque_third_law)
end

function addGeneralizedForcesExternal!(f_generalized::Vector{T}, wrench::Wrench{T}, jac::Nothing,
        torque_third_law::Vector{T}) where {T}

    return nothing
end
function addGeneralizedForcesExternal!(f_generalized::Vector{T}, wrench::Wrench{T}, jac::GeometricJacobian{Matrix{T}},
        torque_third_law::Vector{T}) where {T}

    torque!(torque_third_law, jac, wrench)
    f_generalized .+= torque_third_law
    return nothing
end
