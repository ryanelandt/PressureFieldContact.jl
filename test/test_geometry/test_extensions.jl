
eM_box = eMesh_box()

@testset "extensions" begin
    @test_throws ErrorException volume(as_tri_eMesh(eM_box))
    @test_throws ErrorException area(as_tet_eMesh(eM_box))
end
