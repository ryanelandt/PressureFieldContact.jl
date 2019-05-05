
"""
A structure that is assumed to hold a homogenous transform.
"""
struct basic_dh{T}
    mat::SMatrix{4,4,T,16}
    function basic_dh(t::SVector{3,T}) where {T}
        mat = SMatrix{4,4,T,16}( one(T), zero(T), zero(T), zero(T),
                                zero(T),  one(T), zero(T), zero(T),
                                zero(T), zero(T),  one(T), zero(T),
                                   t[1],    t[2],    t[3],  one(T))
        return new{T}(mat)
    end
    function basic_dh(R::Rotation{3,T}) where {T}
        R = R[:]
        mat = SMatrix{4,4,T,16}(   R[1],    R[2],    R[3], zero(T),
                                   R[4],    R[5],    R[6], zero(T),
                                   R[7],    R[8],    R[9], zero(T),
                                zero(T), zero(T), zero(T),  one(T))
        return new{T}(mat)
    end
    function basic_dh(R::Rotation{3,T}, t::SVector{3,T}) where {T}
        R = R[:]
        mat = SMatrix{4,4,T,16}(   R[1],    R[2],    R[3], zero(T),
                                   R[4],    R[5],    R[6], zero(T),
                                   R[7],    R[8],    R[9], zero(T),
                                   t[1],    t[2],    t[3],  one(T))
        return new{T}(mat)
    end
    function basic_dh(R::SMatrix{3,3,T,9}, t::SVector{3,T}) where {T}
        R = R[:]
        mat = SMatrix{4,4,T,16}(   R[1],    R[2],    R[3], zero(T),
                                   R[4],    R[5],    R[6], zero(T),
                                   R[7],    R[8],    R[9], zero(T),
                                   t[1],    t[2],    t[3],  one(T))
        return new{T}(mat)
    end
    function basic_dh(R::SMatrix{3,3,T,9}) where {T}
        R = R[:]
        mat = SMatrix{4,4,T,16}(   R[1],    R[2],    R[3], zero(T),
                                   R[4],    R[5],    R[6], zero(T),
                                   R[7],    R[8],    R[9], zero(T),
                                zero(T), zero(T), zero(T), one(T))
        return new{T}(mat)
    end
    # function basic_dh(;scale::SVector{3,T}) where {T}
    #     return basic_dh(SMatrix{3,3,T,9}(scale[1], zero(T),  zero(T),
    #                                      zero(T),  scale[2], zero(T),
    #                                      zero(T),  zero(T),  scale[3]))
    # end
    function basic_dh(d::Diagonal{T,SVector{3,T}}) where {T}
        return basic_dh(SMatrix{3,3,T,9}(d))
    end
    basic_dh(s::T) where {T} = basic_dh(s * one(SMatrix{3,3,T,9}))
    basic_dh(mat::SMatrix{4,4,T,16}) where {T} = new{T}(mat)
end

"""
Extracts the rotational and translational part of the transformation matrix.
"""
function dh_R_t(a::basic_dh{T}) where {T}
    a_top = getTop(a.mat)
    t = SVector{3,T}(a_top[10], a_top[11], a_top[12])
    R = SMatrix{3,3,T,9}(a_top[1], a_top[2], a_top[3], a_top[4], a_top[5], a_top[6], a_top[7], a_top[8], a_top[9])
    return R, t
end

function povray_12(a::basic_dh)
    # http://www.povray.org/documentation/view/3.6.0/49/
    # qx = Val00 * px + Val10 * py + Val20 * pz + Val30
    # qy = Val01 * px + Val11 * py + Val21 * pz + Val31
    # qz = Val02 * px + Val12 * py + Val22 * pz + Val32
    #
    # matrix <Val00, Val01, Val02, Val10, Val11, Val12, Val20, Val21, Val22, Val30, Val31, Val32>
    return getTop(a.mat)[:]
end

Base.:*(a::basic_dh, b::basic_dh) = basic_dh(a.mat * b.mat)

function Base.inv(a::basic_dh)
    R, t = dh_R_t(a)
    return basic_dh(R', -R' * t)
end

Base.one(::Type{basic_dh{T}}) where {T} = basic_dh(one(SMatrix{4,4,T,16}))

dh_vector_mul(dh::basic_dh{T1}, v::SVector{3,T2}) where {T1,T2} = unPad(dh.mat * onePad(v))
