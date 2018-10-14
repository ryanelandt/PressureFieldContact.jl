
struct Point4D{V<:AbstractVector}
    v::V
    frame::CartesianFrame3D

    @inline function Point4D(frame::CartesianFrame3D, v::V) where {V}
        @boundscheck length(v) == 4 || throw(DimensionMismatch())
        new{V}(v, frame)
    end
end

struct MatrixTransform{N1,N2,T,N3}
    mat::SMatrix{N1,N2,T,N3}
    from::CartesianFrame3D
    to::CartesianFrame3D
    MatrixTransform(from::CartesianFrame3D, to::CartesianFrame3D, mat::SMatrix{N1,N2,T,N3}) where {N1,N2,T,N3} = new{N1,N2,T,N3}(mat, from, to)
end

Base.:*(t1::Transform3D, t2::MatrixTransform) = begin @framecheck(t1.from, t2.to); MatrixTransform(t2.from, t1.to, mat_mul_SA_bug_circumvent(t1.mat, t2.mat)) end
Base.:*(t1::MatrixTransform, t2::Transform3D) = begin @framecheck(t1.from, t2.to); MatrixTransform(t2.from, t1.to, mat_mul_SA_bug_circumvent(t1.mat, t2.mat)) end
# Base.:*(t1::Transform3D, t2::MatrixTransform) = begin @framecheck(t1.from, t2.to); MatrixTransform(t2.from, t1.to, t1.mat * t2.mat) end
# Base.:*(t1::MatrixTransform, t2::Transform3D) = begin @framecheck(t1.from, t2.to); MatrixTransform(t2.from, t1.to, t1.mat * t2.mat) end
Base.:*(t::MatrixTransform{4,4,T,16}, point::Point3D{SVector{3,T2}}) where {T, T2} = begin @framecheck(t.from, point.frame); Point4D(t.to, t.mat * onePad(point.v)) end
Base.:*(t::MatrixTransform{4,3,T,12}, point::Point3D{SVector{3,T2}}) where {T, T2} = begin @framecheck(t.from, point.frame); Point4D(t.to, t.mat * point.v) end
Base.:*(t::MatrixTransform{3,4,T,12}, point::Point4D{SVector{4,T}}) where {T} = begin @framecheck(t.from, point.frame); Point3D(t.to, t.mat * point.v) end
@inline LinearAlgebra.inv(t::MatrixTransform{N,N,T,N1}) where {N,T,N1} = MatrixTransform(t.to, t.from, inv(t.mat))

getPoint(c::ClippedPolygon{3,T}, k::Int64) where {T} = Point3D(c.frame, c[k])
getPoint(c::ClippedPolygon{4,T}, k::Int64) where {T} = Point4D(c.frame, c[k])

function Base.hcat(p1::Point4D{T}, p2::Point4D{T}, p3::Point4D{T}) where {T}
    @framecheck(p1.frame, p2.frame)
    @framecheck(p2.frame, p3.frame)
    return MatrixTransform(FRAME_ϕ, p1.frame, hcat(p1.v, p2.v, p3.v))
end
