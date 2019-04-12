
@testset "mesh_body_utility" begin
    r = 1.2
    l = 2 * r
    eM_box = output_eMesh_box(r)
    @test area(eM_box) ≈ 6 * l^2
    @test volume(eM_box) ≈ l^3
end
