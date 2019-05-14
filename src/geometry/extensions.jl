
NumericalTricks.volume(eM::eMesh{T1,Nothing}) where {T1} = error("meshes without tets have no volume")
function NumericalTricks.volume(eM::eMesh{T1,Tet}) where {T1}
    vol = 0.0
    for k = 1:n_tet(eM)
        vol += volume(eM.point[eM.tet[k]])
    end
    return vol
end

NumericalTricks.area(eM::eMesh{Nothing,T2}) where {T2} = error("meshes without triangles have no area")
function NumericalTricks.area(eM::eMesh{Tri,T2}) where {T2}
    area_cum = 0.0
    for k = 1:n_tri(eM)
        area_cum += area(vertex_pos_for_tri_ind(eM, k))
    end
    return area_cum
end
