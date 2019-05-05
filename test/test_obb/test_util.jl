
# using StaticArrays
# using Binary_BB_Trees
using .Binary_BB_Trees: findMinSVSV, findMaxSVSV, minMaxToCenterExtent


p1 = SVector(1.0, 5.0, 9.0)
p2 = 2 * p1
p3 = 3 * p1
p4 = 4 * p1

@testset "findMinSVSV/findMaxSVSV" begin
    @test p1 == findMinSVSV(p1, p2) == findMinSVSV(p2, p1)
    @test p1 == findMinSVSV(SVector(p1, p2, p3)) == findMinSVSV(SVector(p3, p2, p1)) == findMinSVSV(SVector(p2, p1, p3))
    @test p1 == findMinSVSV(SVector(p1, p2, p3, p4)) == findMinSVSV(SVector(p4, p3, p2, p1)) == findMinSVSV(SVector(p2, p4, p1, p3))

    @test p4 == findMaxSVSV(p4, p3) == findMaxSVSV(p3, p4)
    @test p4 == findMaxSVSV(SVector(p4, p3, p2)) == findMaxSVSV(SVector(p2, p3, p4)) == findMaxSVSV(SVector(p1, p4, p3))
    @test p4 == findMaxSVSV(SVector(p1, p2, p3, p4)) == findMaxSVSV(SVector(p4, p3, p2, p1)) == findMaxSVSV(SVector(p2, p4, p1, p3))
end

@testset "minMaxToCenterExtent" begin
    center, extent = minMaxToCenterExtent(p1, p2)
    @test center == (p1 + p2) * 0.5
    @test extent == (p2 - p1) * 0.5
end

@testset "calc_min_max" begin
    min_, max_ = calc_min_max(p1, p2)
    @test min_ == p1
    @test max_ == p2
    min_, max_ = calc_min_max(p2, p1)
    @test min_ == p1
    @test max_ == p2
end
