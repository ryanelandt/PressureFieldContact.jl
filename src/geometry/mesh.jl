
abstract type Tri end
abstract type Tet end

"""
$(TYPEDEF)

Data structure for holding geometry data.
"""
struct eMesh{T1<:Union{Nothing,Tri},T2<:Union{Nothing,Tet}}
    point::Vector{SVector{3,Float64}}
    tri::Union{Nothing,Vector{SVector{3,Int64}}}
    tet::Union{Nothing,Vector{SVector{4,Int64}}}
    ϵ::Union{Nothing,Vector{Float64}}
    function eMesh( point::Vector{SVector{3,Float64}},
                    tri::Union{Nothing,Vector{SVector{3,Int64}}},
                    tet::Union{Nothing,Vector{SVector{4,Int64}}}=nothing,
                    ϵ::Union{Nothing,Vector{Float64}}=nothing)

        T1_ = ifelse(tri == nothing, Nothing, Tri)
        T2_ = ifelse(tet == nothing, Nothing, Tet)
        if T2_ == Tet
            @assert(isa(ϵ, Vector{Float64}))
            @assert(length(ϵ) == length(point), "length(ϵ) = $(length(ϵ)) but length(point) = $(length(point))")
            (length(ϵ) != 0) && @assert(0.0 < maximum(ϵ), "normalized penetration extent must be non-negative")
            (length(ϵ) != 0) && @assert(minimum(ϵ) == 0.0, "normalized penetration extent must be zero on the surface of the volume mesh")
            for k = 1:length(tet)
                (0.0 < volume(point[tet[k]])) || error("inverted tetrahedron")
            end
        else
            @assert(ϵ == nothing)
        end
        (T1_ == T2_ == Nothing) && error("a whole lot of nothing")
        return new{T1_,T2_}(point, tri, tet, ϵ)
    end
    function eMesh(hm::HomogenousMesh, tet::Union{Nothing,Vector{SVector{4,Int64}}}=nothing,
            ϵ::Union{Nothing,Vector{Float64}}=nothing)

        point = [SVector{3,Float64}(k) for k = hm.vertices]
        tri = [SVector{3,Int64}(k) for k = hm.faces]
        return eMesh(point, tri, tet, ϵ)
    end
    eMesh{Tri,Nothing}() = eMesh(Vector{SVector{3,Float64}}(), Vector{SVector{3,Int64}}(), nothing, nothing)
    eMesh{Nothing,Tet}() = eMesh(Vector{SVector{3,Float64}}(), nothing, Vector{SVector{4,Int64}}(), Vector{Float64}())
    eMesh{Tri,Tet}()     = eMesh(Vector{SVector{3,Float64}}(), Vector{SVector{3,Int64}}(), Vector{SVector{4,Int64}}(), Vector{Float64}())
end

as_tet_eMesh(e_mesh::eMesh{Tri,Tet}) = eMesh(e_mesh.point, nothing, e_mesh.tet, e_mesh.ϵ)
as_tet_eMesh(e_mesh::eMesh{Nothing,Tet}) = deepcopy(e_mesh)
as_tri_eMesh(e_mesh::eMesh{Tri,Tet}) = eMesh(e_mesh.point, e_mesh.tri, nothing, nothing)
as_tri_eMesh(e_mesh::eMesh{Tri,Nothing}) = deepcopy(e_mesh)

vertex_pos_for_tri_ind(eM::eMesh{Tri,T2}, k::Int64) where {T2} = eM.point[eM.tri[k]]
vertex_pos_for_tet_ind(eM::eMesh{T1,Tet}, k::Int64) where {T1} = eM.point[eM.tet[k]]

n_point(eM::eMesh) = length(eM.point)
n_tri(eM::eMesh) = length(eM.tri)
n_tet(eM::eMesh) = length(eM.tet)
get_tri(eM::eMesh) = eM.tri
get_tet(eM::eMesh) = eM.tet
get_point(eM::eMesh) = eM.point

function Base.empty!(e_mesh::eMesh{T1,T2}) where {T1,T2}
    empty!(e_mesh.point)
    (T1 == Nothing) || empty!(e_mesh.tri)
    (T2 == Nothing) || empty!(e_mesh.tet)
    (T2 == Nothing) || empty!(e_mesh.ϵ)
    return nothing
end

function Base.isempty(e_mesh::eMesh{T1,T2}) where {T1,T2}
    is_emp = isempty(e_mesh.point)
    if T1 == Tri
        is_emp = is_emp || isempty(e_mesh.tri)
    end
    if T2 == Tet
        is_emp = is_emp || isempty(e_mesh.tet)
        is_emp = is_emp || isempty(e_mesh.ϵ)
    end
    return is_emp
end

function Base.append!(eM_1::eMesh{T1,T2}, eM_2::eMesh{T1,T2}) where {T1,T2}
    n_1 = n_point(eM_1)
    if T1 != Nothing
        for k = eM_2.tri
            push!(eM_1.tri, k .+ n_1)
        end
    end
    if T2 != Nothing
        for k = eM_2.tet
            push!(eM_1.tet, k .+ n_1)
        end
        append!(eM_1.ϵ, eM_2.ϵ)
    end
    append!(eM_1.point, eM_2.point)
    return nothing
end

### VERIFICATION

function verify_mesh(eM::eMesh{T1,T2}) where {T1,T2}
    verify_mesh_triangles(eM)
    verify_mesh_tets(eM)
    return nothing
end

verify_mesh_tets(eM::eMesh{T1,Nothing}) where {T1} = nothing
function verify_mesh_tets(eM::eMesh{T1,Tet}) where {T1}
    length_tet = n_tet(eM)
    length_point = n_point(eM)
    length_ϵ = length(eM.ϵ)
    for k = 1:length_tet
        iΔ = eM.tet[k]
        all(1 .<= iΔ .<= length_point) || error("index in tet with sides $(iΔ) not within points")
    end
    (length_ϵ == length_point) || error("number of points not equal to number or ϵ")
    return nothing
end

verify_mesh_triangles(eM::eMesh{Nothing,T2}) where {T2} = nothing
function verify_mesh_triangles(eM::eMesh{Tri,T2}) where {T2}
    length_tri = n_tri(eM)
    length_point = n_point(eM)
    for k = 1:length_tri
        iΔ = eM.tri[k]
        all(1 .<= iΔ .<= length_point) || error("index in triangle with sides $(iΔ) not within points")
    end
    return nothing
end

### MESH MANIPULATION

function eMesh_transform!(e_mesh::eMesh{T1,T2}, extra_args...) where {T1,T2}
    dh = basic_dh(extra_args...)
    point = e_mesh.point
    for k = 1:n_point(e_mesh)
        point[k] = dh_vector_mul(dh, point[k])
    end
end

struct PointObj
    p::SVector{3,Float64}
    w::Float64
    b::Bool
    k::Int64
end

function crop_mesh(eM::eMesh{Tri,Nothing}, plane::SMatrix{1,4,Float64,4})
    function add_cropped_triangle!(eM_new::eMesh{Tri,Nothing}, iΔ::SVector{3,Int64})
        function plane_dot_2(eM_new::eMesh{Tri,Nothing}, plane::SMatrix{1,4,Float64,4}, k::Int64)
            p = get_point(eM_new)[k]
            w = dot(plane, onePad(p))
            b = -1.0e-12 < w
            return PointObj(p, w, b, k)
        end

        o1 = plane_dot_2(eM_new, plane, iΔ[1])
        o2 = plane_dot_2(eM_new, plane, iΔ[2])
        o3 = plane_dot_2(eM_new, plane, iΔ[3])
        s_bool = o1.b + o2.b + o3.b
        (s_bool == 3) && (push!(eM_new.tri, iΔ))  # all entire triangle
        ((s_bool == 0) || (s_bool == 3)) && (return nothing)  # exit function
        (o1.b == o2.b) && ((o1, o2, o3) = (o2, o3, o1))  # put 1 and 2 on opposite sides
        (o2.b == o3.b) && ((o1, o2, o3) = (o3, o1, o2))  # put 2 and 3 on opposite sides
        (o1.b == o2.b) || (o2.b == o3.b) && error("something is wrong")  # sanity check
        push!(eM_new.point, weightPoly(o1.p, o2.p, o1.w, o2.w))
        i12 = length(get_point(eM_new))
        push!(eM_new.point, weightPoly(o2.p, o3.p, o2.w, o3.w))
        i23 = length(get_point(eM_new))
        if o2.b  # add top of triangle
            push!(eM_new.tri, SVector{3,Int64}(i12, o2.k, i23))
        else  # add bottom quadralateral
            push!(eM_new.tri, SVector{3,Int64}(o1.k, i12, i23))
            push!(eM_new.tri, SVector{3,Int64}(o1.k, i23, o3.k))
        end
    end

    eM_new = eMesh{Tri,Nothing}()
    append!(eM_new.point, eM.point)
    for iΔ = eM.tri
        add_cropped_triangle!(eM_new, iΔ)
    end
    mesh_repair!(eM_new)
    return eM_new
end

function invert!(eM::eMesh{Tri,Nothing})
    for k = 1:n_tri(eM)
        tri_k = eM.tri[k]
        eM.tri[k] = SVector{3,Float64}(tri_k[3], tri_k[2], tri_k[1])
    end
end

### MESH REPAIR
function mesh_repair!(e_mesh::eMesh{T1,T2}) where {T1,T2}
    mesh_remove_unused_points!(e_mesh)
    mesh_inplace_rekey!(e_mesh)
    mesh_remove_unused_points!(e_mesh)
    return delete_triangles!(e_mesh)
end

function remove_degenerate!(eM::eMesh)
	area_or_vol(v::SVector{3,SVector{3,Float64}}) = area(v)
	area_or_vol(v::SVector{4,SVector{3,Float64}}) = volume(v)

	tol = 1.0e-6
	for vec_ind = (eM.tet, eM.tri)
		if vec_ind != nothing
			v = [area_or_vol(eM.point[vec_ind[k]]) for k = 1:length(vec_ind)]
			max_v = maximum(v)
			i_bad = findall(v .< max_v * tol)
			deleteat!(vec_ind, i_bad)
		end
	end
end

rekey!(v::Nothing, i::Vector{Int64}) = nothing
rekey!(v::Vector{SVector{N,Int64}}, i::Vector{Int64}) where {N} = replace!(x -> i[x], v)

function mesh_inplace_rekey(e_mesh::eMesh{T1,T2}) where {T1,T2}
    e_mesh_copy = deepcopy(e_mesh)
    mesh_inplace_rekey!(e_mesh_copy)
    return e_mesh_copy
end

function mesh_inplace_rekey!(e_mesh::eMesh{T1,T2}) where {T1,T2}
    # Expresses the mesh using the minumum number of points

    function shortest_side(e_mesh::eMesh)
        shortest_side(v::Nothing) = Inf
        shortest_side(ind::Vector{SVector{N,Int64}}) where {N} =  minimum([shortest_side(e_mesh.point[k]) for k = ind])
        function shortest_side(v::SVector{N,SVector{3,Float64}}) where {N}
            d_min = Inf
            for k = 1:N
                for kk = 1:(k-1)
                    d_min = min(d_min, norm(v[k] - v[kk]))
                end
            end
            return d_min
        end

        return min(shortest_side(e_mesh.tri), shortest_side(e_mesh.tet))
    end

    function inplace_rekey(point::Vector{SVector{3,Float64}}, min_side_length::Float64)
        balltree = BallTree(point; reorder = false)
        idxs = inrange(balltree, point, min_side_length * 0.499)  # Assume that all points closer to each other than hald a side length are duplicates
        return first.(sort!.(idxs))
    end

    if !isempty(e_mesh)  # reducing over an empty collection is not allowed
        min_side_length = shortest_side(e_mesh)
        new_key = inplace_rekey(e_mesh.point, min_side_length)
        rekey!(e_mesh.tri, new_key)
        rekey!(e_mesh.tet, new_key)
        mesh_remove_unused_points!(e_mesh)
        return new_key
    end
    return nothing
end

function mesh_remove_unused_points!(e_mesh::eMesh{T1,T2}) where {T1,T2}
    truth_vector = falses(n_point(e_mesh))
    (T1 != Nothing) && (for k = e_mesh.tri; truth_vector[k] = true; end)
    (T2 != Nothing) && (for k = e_mesh.tet; truth_vector[k] = true; end)
    new_key = cumsum(truth_vector)
    rekey!(e_mesh.tri, new_key)
    rekey!(e_mesh.tet, new_key)
    ind_truth = findall(truth_vector)
    point_new = e_mesh.point[ind_truth]
    empty!(e_mesh.point)
    append!(e_mesh.point, point_new)
    if T2 != Nothing
        ϵ_new = e_mesh.ϵ[ind_truth]
        empty!(e_mesh.ϵ)
        append!(e_mesh.ϵ, ϵ_new)
    end
    return nothing
end

delete_triangles!(e_mesh::eMesh{Nothing,T2}) where {T2} = -9999
function delete_triangles!(e_mesh::eMesh{Tri,T2}) where {T2}
    ### make the first index lowest
    for k = 1:n_tri(e_mesh)
        iΔ = e_mesh.tri[k]
        min_index = minimum(iΔ)
        (min_index == iΔ[1]) || (iΔ = SVector(iΔ[3], iΔ[1], iΔ[2]))
        (min_index == iΔ[1]) || (iΔ = SVector(iΔ[3], iΔ[1], iΔ[2]))
    end

    ### create a dictionary of key repition counts
    key_type = Vector{Tuple{Int64,SVector{3,Int64}}}
    dict_ind = Dict{SVector{3,Int64},key_type}()
    for k = 1:n_tri(e_mesh)
        iΔ = e_mesh.tri[k]
        sort_iΔ = sort(iΔ)
        !haskey(dict_ind, sort_iΔ) && (dict_ind[sort_iΔ] = key_type() )
        push!(dict_ind[sort_iΔ], (k, iΔ))
    end

    ### check for issues is triangle pairs
    i_delete = Vector{Int64}()
    for (key_k, val_k) = dict_ind
        length_val_k = length(val_k)
        if length_val_k == 2
            val_k1 = val_k[1]
            val_k2 = val_k[2]
            (val_k1[2] == val_k2[2]) && error("non-opposing triangles")
            push!(i_delete, val_k1[1], val_k2[1])
        elseif 3 <= length_val_k
            error("something is wrong")
        end
    end

    ### Actually delete the opposing triangles
    i_delete = sort(i_delete)
    deleteat!(e_mesh.tri, i_delete)
	return length(i_delete)
end

function sub_div_mesh(eM_ico::eMesh{Tri,T2}, n_div::Int64) where {T2}
    n_end(n::Int64) = div((n + 1) * n, 2)
    n_start(n::Int64) = 1 + n_end(n - 1)

    function sub_div_triangle(p::SVector{3,SVector{3,Float64}}, n_div::Int64)
        function sub_div_triangle_vert_index(n_div::Int64)
            i_tri = Vector{SVector{3,Int64}}()
            for k = 1:n_div
                for kk = 0:(k - 1)
                    i1 = n_start(k) + kk
                    i2 = i1 + k
                    i3 = i2 + 1
                    push!(i_tri, SVector{3,Int64}(i1, i2, i3))
                end
                for kk = 0:(k - 2)
                    i1 = n_start(k) + kk
                    i2 = i1 + k + 1
                    i3 = i2 - k
                    push!(i_tri, SVector{3,Int64}(i1, i2, i3))
                end
            end
            return i_tri
        end

        function get_new_point(n_vert::Int64, n_div::Int64, p::SVector{3,SVector{3,Float64}})
            function find_layer_1(i_point::Int64)
                i_end_layer = 1
                i_layer = 1
                while i_end_layer < i_point
                    i_layer += 1
                    i_end_layer += i_layer
                end
                norm_extent_1 = ifelse(i_point == 1, 0.0, (i_end_layer - i_point) / (i_layer - 1) )
                return i_layer, norm_extent_1
            end

            i_layer, norm_extent_1 = find_layer_1(n_vert)
            ϕ_1 = (n_div - i_layer + 1) / n_div
            ϕ_2 = (1 - ϕ_1) * norm_extent_1
            ϕ_3 = 1 - ϕ_1 - ϕ_2
            ϕ = SVector{3,Float64}(ϕ_1, ϕ_2, ϕ_3)
            return sum(p .* ϕ)
        end

        point = Vector{SVector{3,Float64}}()
        i_tri_div = sub_div_triangle_vert_index(n_div)
        for k = 1:n_end(n_div + 1)
            push!(point, get_new_point(k, n_div, p))
        end
        return eMesh(point, i_tri_div, nothing, nothing)
    end

    eM_ico_div = eMesh{Tri,Nothing}()
    for k = 1:n_tri(eM_ico)
        p = vertex_pos_for_tri_ind(eM_ico, k)
        append!(eM_ico_div, sub_div_triangle(p, n_div))
    end
    mesh_repair!(eM_ico_div)
    return eM_ico_div
end

### BASIC SHAPES
"""
$(SIGNATURES)

Outputs an `eMesh` for a half-plane.
"""
function output_eMesh_half_plane(plane_w::Float64=1.0, is_include_vis_sides::Bool=false)
    v1, v2, v3 = [SVector{3,Float64}(cos(theta), sin(theta), 0.0) for theta = (0.0, 2*pi/3, 4*pi/3)]
    v4 = SVector{3,Float64}(0.0, 0.0, -1.0) * plane_w
    point = [v1, v2, v3, v4]
    if is_include_vis_sides
        tri = [SVector{3,Int64}(1,2,3), SVector{3,Int64}(1,3,4),  SVector{3,Int64}(1,4,2), SVector{3,Int64}(2,4,3)]
    else
        tri = [SVector{3,Int64}(1,2,3)]
    end
    tet = [SVector{4,Int64}(4,1,2,3)]
    ϵ = [0.0, 0.0, 0.0, plane_w]
    return eMesh(point, tri, tet, ϵ)
end

"""
$(SIGNATURES)

Outputs an `eMesh` for a sphere. Larger values of n_div create finer discretizations.
"""
function output_eMesh_sphere(rad::Union{Float64,SVector{3,Float64}}=1.0, n_div::Int64=4)
    function make_icosahedron()
        φ = Base.MathConstants.golden

        v = Vector{SVector{3,Float64}}()
        for s_1 = (-1.0, +1.0)
            for s_2 = (-1.0, +1.0)
                push!(v, SVector{3,Float64}(    0.0,     s_1, φ * s_2) )
                push!(v, SVector{3,Float64}(    s_1, φ * s_2,     0.0) )
                push!(v, SVector{3,Float64}(φ * s_2,     0.0,     s_1) )
            end
        end

        n_ico_vert = 12
        d = zeros(n_ico_vert, n_ico_vert)
        for k = 1:n_ico_vert
            for kk = 1:n_ico_vert
                d[k, kk] = norm(v[k] - v[kk])
            end
        end

        v_face_vert = Vector{SVector{3,Int64}}()
        b = d .== 2.0
        for i1 = 1:n_ico_vert
            for i2 = 1:n_ico_vert
                for i3 = 1:n_ico_vert
                    if (i1 < i2 < i3) && b[i1, i2] && b[i2, i3] && b[i1, i3]
                        p1 = v[i1]
                        p2 = v[i2]
                        p3 = v[i3]
                        n = triangleNormal(p1, p2, p3)
                        c = normalize(centroid(p1, p2, p3))
                        i_face = ifelse(n ≈ c, SVector{3,Int64}(i1, i2, i3), SVector{3,Int64}(i1, i3, i2))
                        push!(v_face_vert, i_face)
                    end
                end
            end
        end

        push!(v, SVector{3,Float64}(0,0,0))
        v_tet = Vector{SVector{4,Int64}}()
        for k = 1:length(v_face_vert)
            push!(v_tet, SVector{4,Int64}(13, v_face_vert[k]...))
        end
        ϵ = vcat(zeros(n_ico_vert), 1.0)
        return eMesh(v, v_face_vert, v_tet, ϵ)
    end

    function volumize_about(eM::eMesh{Tri,Nothing})
        i_tet = Vector{SVector{4,Int64}}()
        n_vert = n_point(eM)
        n_center = n_vert + 1
        for k = 1:n_tri(eM)
            push!(i_tet, SVector{4,Int64}(n_center, eM.tri[k]...))
        end
        ϵ = vcat(zeros(n_vert), 1.0)
        point = deepcopy(eM.point)
        push!(point, zeros(SVector{3,Float64}))
        return eMesh(point, eM.tri, i_tet, ϵ)
    end

    function project_to_sphere!(eM::eMesh)
        for k = 1:n_point(eM)
            eM.point[k] = normalize(eM.point[k])
        end
        return nothing
    end

    eM_ico = make_icosahedron()
    eM_ico_div = sub_div_mesh(eM_ico, n_div)
    project_to_sphere!(eM_ico_div)

    rad = ones(SVector{3,Float64}) .* rad
    eMesh_transform!(eM_ico_div, Diagonal(rad))

    return volumize_about(eM_ico_div)
end

function output_box_ind()
    oriented_box_faces = (
        SVector{4, Int64}(1,3,5,7),  # -x
        SVector{4, Int64}(2,6,4,8),  # +x
        SVector{4, Int64}(1,5,2,6),  # -y
        SVector{4, Int64}(3,4,7,8),  # +y
        SVector{4, Int64}(1,2,3,4),  # -z
        SVector{4, Int64}(5,7,6,8),  # +z
    )
    ϵ = zeros(Float64,8)
    push!(ϵ, 1.0)
    tri = Vector{SVector{3,Int64}}()
    tet = Vector{SVector{4,Int64}}()
    for k = 1:6
        bf_k = oriented_box_faces[k]
        push!(tri, bf_k[SVector{3,Int64}(1,3,4)])
        push!(tri, bf_k[SVector{3,Int64}(1,4,2)])
    end
    for k = 1:length(tri)
        tri_k = tri[k]
        push!(tet, SVector{4,Int64}(9, tri_k[1], tri_k[2], tri_k[3]))
    end
    return tri, tet, ϵ
end

"""
$(SIGNATURES)

Outputs an `eMesh` for a box.
"""
function output_eMesh_box(r::Union{Float64,SVector{3,Float64}}=1.0, c::SVector{3,Float64}=zeros(SVector{3,Float64}))
    point = [
        SVector{3,Float64}(-1,-1,-1),
        SVector{3,Float64}(+1,-1,-1),
        SVector{3,Float64}(-1,+1,-1),
        SVector{3,Float64}(+1,+1,-1),
        SVector{3,Float64}(-1,-1,+1),
        SVector{3,Float64}(+1,-1,+1),
        SVector{3,Float64}(-1,+1,+1),
        SVector{3,Float64}(+1,+1,+1),
        SVector{3,Float64}( 0, 0, 0),
    ]
    tri, tet, ϵ = output_box_ind()
    e_mesh = eMesh(point, tri, tet, ϵ)
    r = r .* ones(SVector{3,Float64})
    eMesh_transform!(e_mesh, SMatrix{3,3,Float64,9}(r[1], 0, 0, 0, r[2], 0, 0, 0, r[3]))
    eMesh_transform!(e_mesh, c)
    return e_mesh
end

# function make_cone_points(rⁱ::Float64, rᵒ::Float64, h::Float64, k_slice::Int64, tot_slice::Int64)
#     polar_point(r::Float64, θ::Float64, z::Float64) = SVector{3,Float64}(r * cos(θ), r * sin(θ), z)
#
#     @assert(0 < rⁱ < rᵒ)
#     @assert(0 < h)
#     @assert(1 <= k_slice <= tot_slice)
#
#     v = Vector{SVector{3,Float64}}()
#     θ_space = LinRange{Float64}(0.0, 2*pi, tot_slice + 1)
#     θ_1 = θ_space[k_slice]
#     θ_2 = θ_space[k_slice + 1]
#     for k = 1:8
#         θ = ifelse(mod1(k, 4) <= 2, θ_1, θ_2)
#         r = ifelse(isodd(k), rⁱ, rᵒ)
#         z = ifelse(k <= 4, -0.5, +0.5) * h
#         push!(v, polar_point(r, θ, z))
#     end
#     push!(v, polar_point((rⁱ + rᵒ) / 2, (θ_1 + θ_2) / 2, 0.0))
# end
#
# function output_eMesh_slice(rⁱ::Float64, rᵒ::Float64, h::Float64, k_slice::Int64, tot_slice::Int64)
#     point = make_cone_points(rⁱ, rᵒ, h, k_slice, tot_slice)
#     tri, tet, ϵ = output_box_ind()
#     return eMesh(point, tri, tet, ϵ)
# end
#
# function output_eMesh_hole(rⁱ::Float64, rᵒ::Float64, h::Float64, tot_slice::Int64)
#     eM_cone = eMesh{Tri,Tet}()
#     for k = 1:tot_slice
#         append!(eM_cone, output_eMesh_slice(rⁱ, rᵒ, h, k, tot_slice))
#     end
#     return eM_cone
# end


const SF3 = SVector{3,Float64}
const SI3 = SVector{3,Int64}
const SI4 = SVector{4,Int64}

function create_2D_circle(rad::Float64=1.0; n::Int64=12)
    points_2D = Vector{SF3}()
    tri_2D = Vector{SI3}()
    i_center = n + 1
    θ_range = LinRange(0.0, 2 * π, n + 1)[2:end]
    for (k, θ) = enumerate(θ_range)
        s, c = sincos(θ)
        push!(points_2D, rad * SF3(c, s, 0.0))
        push!(tri_2D, SI3(k, mod1(k + 1, n), i_center))
    end
    push!(points_2D, SF3(0, 0, 0))
    return eMesh(points_2D, tri_2D)
end

function smallest_first(v::SI4)
	_, i = findmin(v)
	return SI4(v[i], v[mod1(i + 1, 4)], v[mod1(i + 2, 4)], v[mod1(i + 3, 4)])
end

function add_stuff_for_extrude!(tri::Vector{SI3}, tet::Vector{SI4}, k::Int64, i_tri_k::SI3, n_✴_2D::Int64)
	b1, b2, b3 = i_tri_k
	t4, t5, t6 = i_tri_k + n_✴_2D
	i_center = k + n_✴_2D * 2
	oriented_slice_faces = (
        SI4(b1,b2,t5,t4),
        SI4(b2,b3,t6,t5),
        SI4(b1,t4,t6,b3)
    )
    tri_add = Vector{SI3}()
    push!(tri_add, SI3(b1,b3,b2))
    push!(tri_add, SI3(t4,t5,t6))
    for bf_k = oriented_slice_faces
		bf_k = smallest_first(bf_k)
		push!(tri_add, bf_k[SI3(1,2,3)])
        push!(tri_add, bf_k[SI3(1,3,4)])
    end
    for tri_k = tri_add
        push!(tri, tri_k)
        push!(tet, SI4(i_center, tri_k[1], tri_k[2], tri_k[3]))
    end
end

function verify_trianges_have_same_normal(surf::eMesh{Tri,Nothing})
	n̂ = [triangleNormal(vertex_pos_for_tri_ind(surf, k)) for k = 1:n_tri(surf)]
	n̂_mean = n̂[1]
	for n̂_k = n̂
		(n̂_k ≈ n̂_mean) || error("All triangles must have the same normal.")
	end
	return n̂_mean
end

"""
$(SIGNATURES)

Extrudes a planar `eMesh` of type `eMesh{Tri,Nothing}` into a type `eMesh{Tri,Tet}` in the direction of the plane
normal. This direction is detected automatically. The function errors if the points do not lie on a plane.
"""
function extrude_mesh(surf::eMesh{Tri,Nothing}, thick::Float64)
	n̂ = verify_trianges_have_same_normal(surf)

    tri_2D    = surf.tri
    points_2D = surf.point
	n_▲_2D = length(tri_2D )
	n_✴_2D = length(points_2D)
    point_lo = points_2D .- [n̂ * thick / 2]
    point_hi = points_2D .+ [n̂ * thick / 2]
    point_centroid = [centroid(points_2D[k]) for k = tri_2D]
    point = vcat(point_lo, point_hi, point_centroid)
    ϵ = vcat(zeros(2 * n_✴_2D), ones(n_▲_2D))
    tri = Vector{SI3}()
    tet = Vector{SI4}()
    for (k, i_tri_k) = enumerate(tri_2D)
		add_stuff_for_extrude!(tri, tet, k, i_tri_k, n_✴_2D)
    end
    eM_thick = eMesh(point, tri, tet, ϵ)
	mesh_repair!(eM_thick)
    return eM_thick
end

"""
$(SIGNATURES)

Outputs an `eMesh` for a cylinder.
"""
function output_eMesh_cylinder(rad::Float64=1.0, height::Float64=1.0; n::Int64=6)
	eM_circle = create_2D_circle(rad, n=n)
	return extrude_mesh(eM_circle, height)
end
