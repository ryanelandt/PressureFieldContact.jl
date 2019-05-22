
@testset "tree creation" begin
	eM = eMesh_sphere()

	@test_throws ErrorException eMesh_to_tree(eM)
	num_tri = length(eM.tri)
	num_tet = length(eM.tet)

	tree = eMesh_to_tree(as_tri_eMesh(eM))
	@test leafNumber(tree) == num_tri
	@test Set(extractData(tree)) == Set(1:num_tri)
	@test treeDepth(tree) < (log2(num_tri) * 1.3)

	tree = eMesh_to_tree(as_tet_eMesh(eM))
	@test leafNumber(tree) == num_tet
	@test Set(extractData(tree)) == Set(1:num_tet)
	@test treeDepth(tree) < (log2(num_tet) * 1.3)
end
