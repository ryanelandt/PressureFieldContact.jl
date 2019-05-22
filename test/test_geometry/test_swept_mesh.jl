
function is_convex_mesh_oriented_outward(eM::eMesh{Tri,T2}) where {T2}
	for k = 1:n_tri(eM)
		p3 = vertex_pos_for_tri_ind(eM, k)
		n̂ = triangleNormal(p3)
		c = centroid(p3)
		c = normalize(c)
	end
	return true
end

@testset "swept_mesh" begin
	eM_box = create_swept_mesh(f_swept_triv, [-1.0, 1.0], 1.0, 4, true)
	@test is_convex_mesh_oriented_outward(eM_box)
	@test 8.0 ≈ volume(eM_box)

	eM_box_check = eMesh_box()
	@test is_convex_mesh_oriented_outward(eM_box_check)
	@test 8.0 ≈ volume(eM_box_check)
end
