#
# ### norm²
# function norm²(v::SVector{3,Float64})
#     n2 = muladd(v[1], v[1], 1.0e-66)
#     n2 = muladd(v[2], v[2], n2)
#     return muladd(v[3], v[3], n2)
# end
# function norm²(v::SVector{3,Dual{Type_Tag,Float64,N}}) where {Type_Tag, N}
#     n2 = v[1]^2
#     n2 = muladd(v[2], v[2], n2)
#     return muladd(v[3], v[3], n2)
# end
#
# ### unsafe_norm
# @inline unsafe_norm(v::SVector{3,T}) where {T} = sqrt(norm²(v))
#
# ### safe_norm
# @inline safe_norm(v::SVector{3,Float64}) = unsafe_norm(v)
# function safe_norm(v::SVector{3,Dual{Type_Tag,Float64,N}}) where {Type_Tag, N}
#     n = unsafe_norm(v)
#     return ifelse(value(n) == 0.0, zero(Dual{Type_Tag,Float64,N}), n)
# end
#
# ### unsafe_inv_norm
# @inline unsafe_inv_norm(v::SVector{3,T}) where {T} = 1.0 / unsafe_norm(v)
#
# ### safe_inv_norm
# function safe_inv_norm(v::SVector{3,Float64})
#     n = unsafe_norm(v)
#     return ifelse(n == 0.0, 0.0, 1.0 / n)
# end
# function safe_inv_norm(v::SVector{3,Dual{Type_Tag,Float64,N}}) where {Type_Tag, N}
#     n = unsafe_norm(v)
#     return ifelse(value(n) == 0.0, zero(Dual{Type_Tag,Float64,N}), 1.0 / n)
# end
#
# ### safe_normalize
# safe_normalize(v::SVector{3,T}) where {T} = v * safe_inv_norm(v)
#
# ### unsafe_normalize
# unsafe_normalize(v::SVector{3,T}) where {T} = v * unsafe_inv_norm(v)
#
# ### safe_inv_norm²
# function safe_inv_norm²(v::SVector{3,T}) where {T}
#     n2 = norm²(v)
#     return ifelse(n2 == 0.0, zero(T), 1.0 / n2)
# end
#
# ### safe_scalar_divide
# function safe_scalar_divide(a::T, b::T) where {T}
#     (b != 0.0) && (return a / b)
#     (a == 0.0) && (return zero(T))
#     error("attempted to divide non-zero $a by zero $b you should not do this")
# end
