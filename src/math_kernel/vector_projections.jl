vec_proj(vec_in::SVector{3, TD}, n̂::SVector{3,TD}) where {TD} = n̂ * dot(n̂, vec_in)
function vec_sub_vec_proj(v::SVector{3,T}, n̂::SVector{3,T}) where T
	t = -dot(v, n̂)
	return SVector{3,T}(muladd(t, n̂[1], v[1]),
						muladd(t, n̂[2], v[2]),
						muladd(t, n̂[3], v[3]))
end

function a_dot_one_pad_b(a::Union{SMatrix{1,4,T1,4},SVector{4,T1}}, b::SVector{3,T2}) where {T1,T2}
	d =    muladd(a[1], b[1], a[4])
	d =    muladd(a[2], b[2], d)
	return muladd(a[3], b[3], d)
end
