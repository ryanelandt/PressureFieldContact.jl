
"""
Stores a convex planar polygon in with at most 8 sides. The polygon should be in 3 dimensions (Cartesian coordinates)
or 4 dimensions (tetrahedral coondinates).
"""
struct poly_eight{N,T}
    n::Int64
    v::NTuple{8,SVector{N,T}}
    function poly_eight{N,T}() where {N,T}
        return new{N,T}(0)
    end
    function poly_eight(n::Int64, v::NTuple{8,SVector{N,T}}) where {N,T}
        return new{N,T}(n, v)
    end
    function poly_eight(v::NTuple{3,SVector{N,T}}) where {N,T}
        return new{N,T}(3, (v[1], v[2], v[3], v[1], v[1], v[1], v[1], v[1]))
    end
    function poly_eight(v::NTuple{4,SVector{N,T}}) where {N,T}
        return new{N,T}(4, (v[1], v[2], v[3], v[4], v[1], v[1], v[1], v[1]))
    end
end

@inline Base.isempty(p::poly_eight) = (p.n == 0)
@inline Base.length(p::poly_eight) = p.n
@inline Base.getindex(p::poly_eight, k::Int64) = p.v[k]

"""
Finds the centroid of a 3 dimensional poly_eight by dividing the area into triangles. This funciton requires a reference
normal n̂.
"""
function NumericalTricks.centroid(p_new::poly_eight{3,T}, n̂::SVector{3,T}) where {T}
	cart_a = p_new.v[1]
	cart_c = p_new.v[2]  # because c becomes b
	cum_sum  = zero(T)
	cum_prod = zeros(SVector{3,T})
	for k = 3:p_new.n
		cart_b = cart_c
		cart_c = p_new.v[k]
		area_ = triangle_area((cart_a, cart_b, cart_c), n̂)
		cum_prod += area_ * centroid(cart_a, cart_b, cart_c)
		cum_sum  += area_
	end
	if cum_sum == 0.0  # just return an arbitrary vertex if the area is zero
		return cum_sum, cart_a
	else
		return cum_sum, cum_prod / cum_sum
	end
end

"""
One pads each element of a 3 dimensional polygon and then multiplies it by a 4x4 matrix m. This funciton is used to
convert each element of a Cartesian polygon into tetrahedral coordinates.
"""
function NumericalTricks.one_pad_then_mul(m::SMatrix{4,4,T1,16}, p::poly_eight{3,T2}) where {T1,T2}
    t1 = one_pad_then_mul(m, p[1])
    t2 = one_pad_then_mul(m, p[2])
    t3 = one_pad_then_mul(m, p[3])
    t4 = one_pad_then_mul(m, p[4])
    length_p = length(p)
    if length_p <= 4
        return poly_eight(length_p, (t1, t2, t3, t4, t1, t1, t1, t1))
    else
        t5 = one_pad_then_mul(m, p[5])
        t6 = one_pad_then_mul(m, p[6])
        t7 = one_pad_then_mul(m, p[7])
        t8 = one_pad_then_mul(m, p[8])
        return poly_eight(length_p, (t1, t2, t3, t4, t5, t6, t7, t8))
    end
end

"""
Multiplies a 4 dimensional polygon by a 4x4 matrix and then unpads the 4 dimensional polygon. This function is used to
convert a polygon in tetrahedral coordinates to a polygon in Cartesian coordinates.
"""
function NumericalTricks.mul_then_un_pad(m::SMatrix{4,4,T1,16}, p::poly_eight{4,T2}) where {T1,T2}
    t1 = mul_then_un_pad(m, p[1])
    t2 = mul_then_un_pad(m, p[2])
    t3 = mul_then_un_pad(m, p[3])
    t4 = mul_then_un_pad(m, p[4])
    length_p = length(p)
    if length_p <= 4
        return poly_eight(length_p, (t1, t2, t3, t4, t1, t1, t1, t1))
    else
        t5 = mul_then_un_pad(m, p[5])
        t6 = mul_then_un_pad(m, p[6])
        t7 = mul_then_un_pad(m, p[7])
        t8 = mul_then_un_pad(m, p[8])
        return poly_eight(length_p, (t1, t2, t3, t4, t5, t6, t7, t8))
    end
end

"""
Sets small coondinates of a polygon in tetrahedral coordinates to zero. This function addresses a specific probability
one degeneracy that results when clipping an edge in the same place twice.
"""
function zero_small_coordinates(p::poly_eight{4,T}) where {T}
    function zero_if_small(x::SVector{4,T}) where {T}
        abs_x = 1.0e-14 .< abs.(value.(x))
        return x .* abs_x
    end

    v1 = zero_if_small(p[1])
    v2 = zero_if_small(p[2])
    v3 = zero_if_small(p[3])
    v4 = zero_if_small(p[4])
    length_p = length(p)
    if length_p <= 4
        return poly_eight(length_p, (v1, v2, v3, v4, v1, v1, v1, v1))
    else
        v5 = zero_if_small(p[5])
        v6 = zero_if_small(p[6])
        v7 = zero_if_small(p[7])
        v8 = zero_if_small(p[8])
        return return poly_eight(length_p, (v1, v2, v3, v4, v5, v6, v7, v8))
    end
end
