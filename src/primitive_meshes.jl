
function create_surface_half_plane(plane_w::Float64=1.0)
    # R_z = RotZ(2*pi/3)
    # v1 = SVector{3,Float64}(1.0, 0.0, 0.0) * plane_w
    # v2 = R_z * v1
    # v3 = R_z * v2
    v1, v2, v3 = [SVector{3,Float64}(cos(theta), sin(theta), 0.0) for theta = (0.0, 2*pi/3, 4*pi/3)]
    v4 = SVector{3,Float64}(0.0, 0.0, -1.0) * plane_w
    plane_vertices = [v1, v2, v3, v4]
    plane_tri_ind = [SVector{3,Int64}(1,2,3)]
    return HomogenousMesh(vertices=plane_vertices, faces=plane_tri_ind)
end

function create_volume_half_plane(plane_w::Float64=1.0; Ē::Float64=1.0e6, d::Float64=0.1, μ::Float64=0.3, χ::Float64=1.0)
    plane_c_prop = ContactProperties(Ē=Ē, μ=μ, χ=χ, d=d)
    plane_h_mesh = create_surface_half_plane(plane_w)
    plane_tet_ind = [SVector{4,Int64}(1,2,4,3)]
    plane_ϵ = [0.0, 0.0, 0.0, -plane_w * plane_c_prop.d⁻¹]
    plane_tet_mesh = TetMesh(plane_h_mesh.vertices, plane_tet_ind, plane_ϵ, plane_c_prop)
    return plane_h_mesh, plane_tet_mesh
end

function create_surface_box(rad::Union{Float64,SVector{3,Float64}}=1.0;
        center::SVector{3,Float64}=SVector{3,Float64}(0,0,0))

    rad = rad .* SVector{3,Float64}(1,1,1)  # output is always SVector{3,Float64}
    return HomogenousMesh(vertices=sizeCenterBoxPoints(rad, center), faces=outputBoxTriInd())
end

function create_volume_box(rad::Union{Float64,SVector{3,Float64}}=1.0;
        center::SVector{3,Float64}=SVector{3,Float64}(0,0,0), fill_frac::Float64=1.0, Ē::Float64=1.0e6,
        μ::Float64=0.3, χ::Float64=1.0)

    rad = rad .* SVector{3,Float64}(1,1,1)  # output is always SVector{3,Float64}
    (0.01 <= fill_frac <= 1.0) || error("foam depth needs to be between 0.01 and 1.0 of minimum box radius")
    d = minimum(rad) * fill_frac
    rad_inner = rad - d

    points_16 = vcat(sizeCenterBoxPoints(rad, center), sizeCenterBoxPoints(rad_inner, center))
    dict_points = determineUniquePoints(points_16)
    vec_nondegenerate_tet = Vector{SVector{4,Int64}}()
    for ind_k = outputDividedCube()
        ss_k = SVector{4,Int64}([dict_points[k] for k = points_16[ind_k]])
        (length(unique(ss_k)) == 4) && push!(vec_nondegenerate_tet, ss_k)
    end
    points_12 = Vector{SVector{3,Float64}}(undef, length(dict_points))
    for (key_k, k) = dict_points
        points_12[k] = key_k
    end
    ϵ = [ifelse(k <= 8, 0.0, -1.0) for k = 1:length(points_12)]

    box_h_mesh = HomogenousMesh(vertices=points_12, faces=outputBoxTriInd())
    box_tet_mesh = TetMesh(points_12, vec_nondegenerate_tet, ϵ, ContactProperties(Ē=Ē, μ=μ, χ=χ, d=d))
    return box_h_mesh, box_tet_mesh
end

is_good_orientation(v) = (v[1] == 1) || (v[1] == 8)
rotateIndexAboutZ(v) = v[[3, 1, 4, 2, 7, 5, 8, 6]]
cycleArray!(v) = pushfirst!(v, pop!(v))

function makeMinFirst!(ix)
    cycleArray!(ix)
    (findmin(ix)[2] == 1) && (return nothing)
    makeMinFirst!(ix)
end

function outputOrientedBoxFaces()
    return [SVector{4, Int64}(1,3,5,7),  # -x
            SVector{4, Int64}(2,6,4,8),  # +x
            SVector{4, Int64}(1,5,2,6),  # -y
            SVector{4, Int64}(3,4,7,8),  # +y
            SVector{4, Int64}(1,2,3,4),  # -z
            SVector{4, Int64}(5,7,6,8)]  # +z
end

function outputDividedCube()
    a6 = Vector{SVector{4,Int64}}()
    i_face_all = outputOrientedBoxFaces()
    for k_cube_side = 1:6
        i_face_k = Vector{Int64}(i_face_all[k_cube_side])
        v = vcat(i_face_k, i_face_k .+ 8)
        !is_good_orientation(v) && (v = rotateIndexAboutZ(v))
        !is_good_orientation(v) && (v = rotateIndexAboutZ(v))
        !is_good_orientation(v) && (v = rotateIndexAboutZ(v))
        !is_good_orientation(v) && error("should be centered now")
        magic_num = v[1]
        for s_xyz = ([2, 4, 8, 6], [4, 3, 7, 8], [5, 6, 8, 7])
            ix = v[s_xyz]
            makeMinFirst!(ix)
            push!(a6, SVector{4,Int64}(magic_num, ix[1], ix[2], ix[3]))
            push!(a6, SVector{4,Int64}(magic_num, ix[1], ix[3], ix[4]))
        end
    end
    return a6
end

function determineUniquePoints(points_16)
    i_unique = 0
    dict_points = Dict{SVector{3,Float64},Int64}()
    for (k, point_k) = enumerate(points_16)
        if !haskey(dict_points, point_k)
            i_unique += 1
            dict_points[point_k] = i_unique
        end
    end
    return dict_points
end

function outputBoxTriInd()
    v_tri = Vector{SVector{3,Int64}}()
    [twoTriangles!(v_tri, k) for k = outputOrientedBoxFaces()]
    return v_tri
end

box_indices() = [SVector(Tuple(k_cart)) for k_cart = CartesianIndices((2,2,2))[:]]
basicBoxPoints() = [SVector{3,Float64}(ifelse.(Tuple(k).==1,-1.0,+1.0)) for k = box_indices()]

function sizeCenterBoxPoints(scale::SVector{3,Float64}=SVector{3,Float64}(1,1,1), center::SVector{3,Float64}=SVector{3,Float64}(0,0,0))
    v = basicBoxPoints()
    for k = eachindex(v)
        v[k] = v[k] .* scale .+ center
    end
    return v
end

function twoTriangles!(v_tri, i_in)
    push!(v_tri, i_in[SVector{3, Int64}(1,3,4)])
    push!(v_tri, i_in[SVector{3, Int64}(1,4,2)])
    return nothing
end

# function outputBoxVolMesh(;center::SVector{3,Float64}=SVector{3,Float64}(0,0,0),
#                           fill_frac::SVector{3,Float64}=SVector{3,Float64}(1,1,1),
#                           rad::SVector{3,Float64}=SVector{3,Float64}(1,1,1))
#
#     all(0.01 .<= fill_frac .<= 1.0) || error("foam depth needs to be between 0.01 and 1.0 of box radius")
#     rad_inner = rad .* (1 .- fill_frac)
#     points_16 = vcat(sizeCenterBoxPoints(rad, center), sizeCenterBoxPoints(rad_inner, center))
#     dict_points = determineUniquePoints(points_16)
#     vec_nondegenerate_tet = Vector{SVector{4,Int64}}()
#     for ind_k = outputDividedCube()
#         ss_k = SVector{4,Int64}([dict_points[k] for k = points_16[ind_k]])
#         (length(unique(ss_k)) == 4) && push!(vec_nondegenerate_tet, ss_k)
#     end
#     points_12 = Vector{SVector{3,Float64}}(undef, length(dict_points))
#     for (key_k, k) = dict_points
#         points_12[k] = key_k
#     end
#     strain = [ifelse(k <= 8, 0.0, -1.0) for k = 1:length(points_12)]
#
#     mesh_box = HomogenousMesh(vertices=points_12, faces=outputBoxTriInd())
#     return mesh_box, vec_nondegenerate_tet, strain
# end
