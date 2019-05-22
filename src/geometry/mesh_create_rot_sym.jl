
function add_triangle!(point::Vector{SVector{3,Float64}}, i::Vector{SVector{3,Int64}}, a::SVector{3,Float64},
        b::SVector{3,Float64}, c::SVector{3,Float64})

    push!(i, SVector{3,Int64}(1, 2, 3) .+ length(point))
    push!(point, a, b, c)
    return nothing
end

function add_rectangle!(point::Vector{SVector{3,Float64}}, i::Vector{SVector{3,Int64}}, a⁻::SVector{3,Float64},
        a⁺::SVector{3,Float64}, b⁻::SVector{3,Float64}, b⁺::SVector{3,Float64})

    add_triangle!(point, i, a⁻, a⁺, b⁺)
    add_triangle!(point, i, a⁻, b⁺, b⁻)
end

function add_slice(point::Vector{SVector{3,Float64}}, i::Vector{SVector{3,Int64}}, v::NTuple{N,SVector{3,Float64}},
        θ⁰::Float64, θ¹::Float64) where {N}
    # adds a rotationally symmetric slice of point sequence v to the vector of points (point) and vector of indices i

    function add_rot!(point::Vector{SVector{3,Float64}}, i::Vector{SVector{3,Int64}}, a::SVector{3,Float64},
            b::SVector{3,Float64}, RZ_0::RotZ{Float64}, RZ_1::RotZ{Float64})

        a⁻ = RZ_0 * a
        a⁺ = RZ_1 * a
        b⁻ = RZ_0 * b
        b⁺ = RZ_1 * b

        # TODO: is this necessary given mesh repair_mesh?
        if (a[1] == 0.0)  # && (a[2] == 0.0)  # point has no xy extent
            add_triangle!(point, i, a, b⁺, b⁻)
        elseif (b[1] == 0.0)  # && (b[2] == 0.0)
            add_triangle!(point, i, a⁻, a⁺, b)
        else
            add_rectangle!(point, i, a⁻, a⁺, b⁻, b⁺)
        end
    end

    RZ_0 = RotZ(θ⁰)
    RZ_1 = RotZ(θ¹)
    for k = 1:(N - 1)
        add_rot!(point, i, v[k], v[k + 1], RZ_0, RZ_1)
    end
end

function obj_from_point_sequence(point_vec_2D::Vector{SVector{2,Float64}}, n_theta::Int64=10)
    # creates a rotationally symmetric object by rotating a vector of 2D points (x,z) about the z axis to create a
    # vector of 3D points (x,y,z) 

    ### FIXES degenerate points near the axis of rotation
    for k = 1:length(point_vec_2D)
        item_1, item_2 = point_vec_2D[k]
        tol = 1.0e-12
        if item_1 <= -tol
            error("negative extent")
        elseif 0.0 < item_1 <= tol
            @warn("point near rotation axis set to 0.0")
            point_vec_2D[k] = SVector{2,Float64}(0.0, item_2)
        end
    end

    point = Vector{SVector{3,Float64}}()
    i = Vector{SVector{3,Int64}}()
    point_vec = [SVector{3,Float64}(k[1], 0.0, k[2]) for k = point_vec_2D]

    theta = collect(LinRange{Float64}(0, 2*pi, n_theta + 1)) .+ (pi/2)
    for k = 1:n_theta
        add_slice(point, i, Tuple(point_vec), theta[k], theta[mod1(k + 1, n_theta)])
    end

    hex_mesh = eMesh(point, i, nothing, nothing)
    mesh_repair!(hex_mesh)
    return hex_mesh
end
