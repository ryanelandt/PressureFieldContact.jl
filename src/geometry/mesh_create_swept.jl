
"""
$(SIGNATURES)

Circular path for a swept mesh.
Specify the radius by evaluating it ahead of time.
For example, define `my_circle(θ) = f_swept_circle(0.1, θ)`, and then pass `my_circle` to [`create_swept_mesh`](@ref).
"""
function f_swept_circle(r::Float64, θ::Float64)
    n̂_1 = SVector{3,Float64}( cos(θ), sin(θ), 0.0)
    n̂_2 = SVector{3,Float64}(-sin(θ), cos(θ), 0.0)
    return r * n̂_1, n̂_1, n̂_2
end

"""
$(SIGNATURES)

Straight path for a swept mesh.
Pass this to [`create_swept_mesh`](@ref) to make a basic swept mesh.
"""
function f_swept_triv(θ::Float64)
	n̂_1 = SVector(0.0, 0.0, -1.0)
	n̂_2 = SVector(0.0, 1.0,  0.0)
	return n̂_2 * θ, n̂_1, n̂_2
end

function add_rot_sym_segment!(eM, fun_gen, θ, ϕ, rad, is_open::NTuple{2,Bool})
    # x̂ is in the "radial" direction
    # ŷ is along the path

    p1, x̂1, ŷ1 = fun_gen(θ[1])
    p2, x̂2, ŷ2 = fun_gen(θ[2])
    p3 = (p1 + p2) * 0.5
    p4 = p1 + AngleAxis(ϕ[1], ŷ1...) * x̂1 * rad[1]
    p6 = p1 + AngleAxis(ϕ[2], ŷ1...) * x̂1 * rad[1]
    p5 = p2 + AngleAxis(ϕ[1], ŷ2...) * x̂2 * rad[2]
    p7 = p2 + AngleAxis(ϕ[2], ŷ2...) * x̂2 * rad[2]
    n_offset = length(eM.point)
    i1 = SVector{4,Int64}(1,3,4,6) + n_offset
    i2 = SVector{4,Int64}(3,2,5,7) + n_offset
    i3 = SVector{4,Int64}(3,4,6,7) + n_offset
    i4 = SVector{4,Int64}(4,3,5,7) + n_offset
    i5 = SVector{3,Int64}(4,6,7) + n_offset
    i6 = SVector{3,Int64}(4,7,5) + n_offset
    append!(eM.point, [p1, p2, p3, p4, p5, p6, p7])
    append!(eM.tri, [i5, i6])
    append!(eM.tet, [i1, i2, i3, i4])
    ϵ = zeros(7)
    ϵ[1:3] .= 1.0
    if is_open[1]
        ϵ[1] = 0.0
        push!(eM.tri, SVector{3,Int64}(1,6,4) + n_offset)
    end
    if is_open[2]
        ϵ[2] = 0.0
        push!(eM.tri, SVector{3,Int64}(2,5,7) + n_offset)
    end
    append!(eM.ϵ, ϵ)
    return nothing
end

"""
$(SIGNATURES)

Creates a mesh by sweeping a 3D path for the given inputs:
* `fun_gen`: function that takes a single input arc length and outputs 1.) the position on path, 2.) a single direction normal to path and 3.) direction along the path.
* `lr`: arc length locations of nodes on path
* `rad`: thickness of swept mesh
* `n_side`: number of sides of swept path
* `is_open`: if the path starts and ends at the same point (e.g. ring) set this to false
* `rot_half`: if the sides of the path are "off" set this to false
"""
function create_swept_mesh(fun_gen::Function, lr::Union{LinRange,Vector{Float64}}, rad::Float64, n_side::Int64=4, is_open::Bool=true; rot_half::Bool=true)
	l_lr = length(lr)
	if isa(rad, Float64)
		rad = zeros(l_lr) .+ rad
	else
		(l_lr == length(rad)) || error("the length of lr and length of rad must be the same")
	end
    eM = eMesh{Tri,Tet}()
    Δ_ϕ = 2 * pi / n_side
    rad = rad ./ cos(Δ_ϕ / 2)
	n_θ = l_lr - 1
    for k_θ = 1:n_θ
        θ = (lr[k_θ], lr[k_θ + 1])
		rad_k = (rad[k_θ], rad[k_θ + 1])
        for k_ϕ = 1:n_side
            ϕ_0 = Δ_ϕ * (k_ϕ  - 0.5 * rot_half)
            ϕ_1 = ϕ_0 + Δ_ϕ
            ϕ = (ϕ_0, ϕ_1)
            if is_open
                is_open_ = (k_θ==1, k_θ==n_θ)
            else
                is_open_ = (false, false)
            end
            add_rot_sym_segment!(eM, fun_gen, θ, ϕ, rad_k, is_open_)
        end
    end
	remove_degenerate!(eM)
    mesh_repair!(eM)
    return eM
end

function f_swept_helix(θ::Float64, coil_sep::Float64)
    delta_z = 1 / (2 * pi) * coil_sep
    r = SVector{3,Float64}(cos(θ), sin(θ), θ * delta_z)
    dir_1 = normalize(SVector{3,Float64}(-sin(θ), cos(θ), delta_z))
    dir_2 = SVector{3,Float64}( cos(θ), sin(θ), 0.0)
    return r, dir_1, dir_2
end

f_swept_circle(θ::Float64) = f_swept_helix(θ, 0.0)
