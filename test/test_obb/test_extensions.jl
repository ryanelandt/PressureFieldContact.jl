
@testset "extensions" begin
    c = SVector(99.0, 99.0, 99.0)
    e = SVector(1.0, 2.0, 3.0)

    aabb = OBB(c, e, one(SMatrix{3,3,Float64,9}))
    @test 8 * 6.0 == volume(aabb)
    @test 4 * 22.0 == area(aabb)
end
