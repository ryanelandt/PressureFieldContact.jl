
p1 = SVector{2,Float64}(0.0, 0.0)
p2 = SVector{2,Float64}(1.0, 0.0)
p3 = SVector{2,Float64}(1.0, 1.0)
p4 = SVector{2,Float64}(0.0, 1.0)
v = [p1, p2, p3, p4]


@testset "rotationally symmetric surface mesh creation" begin
    mesh_out = obj_from_point_sequence(v, 4)
    @test length(mesh_out.point) == 10
end
