
v1 = SVector{3,Float64}(0.0, 0.0, 0.0)
v2 = SVector{3,Float64}(1.0, 0.0, 0.0)
v3 = SVector{3,Float64}(0.0, 1.0, 0.0)
v4 = SVector{3,Float64}(0.0, 0.0, 1.0)

@testset "triangle" begin
    tri_Sep = (v1, v2, v3)
    tri_SV  = (SVector{3,SVector{3,Float64}}(tri_Sep...),)
    tri_Tup  = ((tri_Sep...),)
    for tri_rep = (tri_Sep, tri_SV, tri_Tup)
        @test area(tri_rep...) ≈ 0.5
        @test centroid(tri_rep...) ≈ SVector{3,Float64}(1/3,1/3,0)
        @test triangleNormal(tri_rep...) ≈ SVector{3,Float64}(0,0,1)
    end
    @test triangle_area(tri_Sep, triangleNormal(tri_Sep...)) ≈ 0.5
end

@testset "tetrahedron" begin
    tet_Tup = (v1, v2, v3, v4)
    tet_SV  = (SVector{4,SVector{3,Float64}}(tet_Tup...),)
    for tet_rep = (tet_SV, tet_Tup)
        @test volume(tet_rep...) ≈ 1/6
        @test centroid(tet_rep...) ≈ SVector{3,Float64}(1,1,1)/4
    end
end
