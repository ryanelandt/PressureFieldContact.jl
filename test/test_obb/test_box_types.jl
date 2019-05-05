

# @testset "AABB" begin
#     c = SVector(99.0, 99.0, 99.0)
#     e = SVector(1.0, 2.0, 3.0)
#
#     aabb = AABB(c, e)
#
#     @test aabb isa AABB
#     @test aabb.c == c
#     @test aabb.e == e
#     @test AABB(aabb) == aabb
# end

@testset "OBB" begin
    c = SVector(99.0, 99.0, 99.0)
    e = SVector(1.0, 2.0, 3.0)
    R = rand(RotMatrix{3})[:,:]

    obb = OBB(c, e, R)
    # aabb = AABB(c, e)
    # obb_via_aabb = OBB(aabb)

    @test obb isa OBB
    @test obb.c == c
    @test obb.e == e
    @test obb.R == R
    # @test obb_via_aabb.R == I
    # @test aabb == AABB(aabb, aabb)
    # @test obb_via_aabb == OBB(obb_via_aabb, obb_via_aabb)
end
