unit_cube_points = SoftContact.basicBoxPoints()

@testset "unit cube" begin
    @test unit_cube_points[1] == SVector{3,Float64}(-1.0, -1.0, -1.0)
    @test unit_cube_points[2] == SVector{3,Float64}(+1.0, -1.0, -1.0)
    @test unit_cube_points[3] == SVector{3,Float64}(-1.0, +1.0, -1.0)
    @test unit_cube_points[4] == SVector{3,Float64}(+1.0, +1.0, -1.0)
    @test unit_cube_points[5] == SVector{3,Float64}(-1.0, -1.0, +1.0)
    @test unit_cube_points[6] == SVector{3,Float64}(+1.0, -1.0, +1.0)
    @test unit_cube_points[7] == SVector{3,Float64}(-1.0, +1.0, +1.0)
    @test unit_cube_points[8] == SVector{3,Float64}(+1.0, +1.0, +1.0)
end

@testset "orientation: cube face / triangle" begin
    i_box_faces = SoftContact.outputOrientedBoxFaces()
    tup_dir = (
        SVector{3,Float64}(-1.0,  0.0,  0.0),
        SVector{3,Float64}(+1.0,  0.0,  0.0),
        SVector{3,Float64}( 0.0, -1.0,  0.0),
        SVector{3,Float64}( 0.0, +1.0,  0.0),
        SVector{3,Float64}( 0.0,  0.0, -1.0),
        SVector{3,Float64}( 0.0,  0.0, +1.0)
    )
    for (k, dir_k) = enumerate(tup_dir)
        v_tri = Vector{SVector{3,SVector{3,Float64}}}()
        SoftContact.twoTriangles!(v_tri, unit_cube_points[i_box_faces[k]])
        tri_1 = v_tri[1]
        tri_2 = v_tri[2]
        @test dir_k ≈ triangleNormal(tri_2)
        @test dir_k ≈ triangleNormal(tri_1)
        area_tri_1 = area(tri_1)
        area_tri_2 = area(tri_2)
        @test area_tri_1 ≈ area_tri_2
        @test 4.0 ≈ (area_tri_1 + area_tri_2)
        centroid_1 = centroid(tri_1)
        centroid_2 = centroid(tri_2)
        @test dir_k ≈ ((centroid_1 + centroid_2) / 2)
    end
end
