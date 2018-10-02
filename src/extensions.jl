NumericalTricks.safe_normalize(p::FreeVector3D{SVector{3,T}}) where {T} = FreeVector3D(p.frame, safe_normalize(p.v))

Tri_Tet_Intersections.area(p1::Point3D{T}, p2::Point3D{T}, p3::Point3D{T}) where {T} = area(p1.v, p2.v, p3.v)
Tri_Tet_Intersections.getTop(m::MatrixTransform{4,N2,T,N3}) where {N2,T,N3} = MatrixTransform(m.from, m.to, getTop(m.mat))
function Tri_Tet_Intersections.add!(c::ClippedPolygon{4,T}, p::Point4D{SVector{4,T}}) where {T}
    @framecheck(c.frame, p.frame)
    add!(c, p.v)
    return nothing
end

function Tri_Tet_Intersections.tet_clip_poly_to_cartesian!(poly_3D::ClippedPolygon{3,T}, poly_4D_1::ClippedPolygon{4,T}, A_w_zeta_top::MatrixTransform) where{T}
    @framecheck(poly_4D_1.frame, A_w_zeta_top.from)
    @framecheck(poly_3D.frame, A_w_zeta_top.to)
    tet_clip_poly_to_cartesian!(poly_3D, poly_4D_1, A_w_zeta_top.mat)
    return nothing
end
