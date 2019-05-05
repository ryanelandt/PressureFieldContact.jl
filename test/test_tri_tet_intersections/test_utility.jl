
@testset "utility" begin
    p1 = rand(SVector{3,Float64})
    p2 = rand(SVector{3,Float64})
    p3 = rand(SVector{3,Float64})
    p4 = rand(SVector{3,Float64})

    m3 = asMat(p1, p2, p3)
    m4 = asMat(p1, p2, p3, p4)
    m3_pad = asMatOnePad(p1, p2, p3)
    m4_pad = asMatOnePad(p1, p2, p3, p4)
    @test all(m3 .== hcat(p1, p2, p3))
    @test all(m4 .== hcat(p1, p2, p3, p4))
    @test all(m3_pad .== hcat(onePad(p1), onePad(p2), onePad(p3)))
    @test all(m4_pad .== hcat(onePad(p1), onePad(p2), onePad(p3), onePad(p4)))

    m3 = asMat( (p1, p2, p3) )
    m4 = asMat( (p1, p2, p3, p4) )
    m3_pad = asMatOnePad( (p1, p2, p3) )
    m4_pad = asMatOnePad( (p1, p2, p3, p4) )
    @test all(m3 .== hcat(p1, p2, p3))
    @test all(m4 .== hcat(p1, p2, p3, p4))
    @test all(m3_pad .== hcat(onePad(p1), onePad(p2), onePad(p3)))
    @test all(m4_pad .== hcat(onePad(p1), onePad(p2), onePad(p3), onePad(p4)))

    m3 = asMat( SVector(p1, p2, p3) )
    m4 = asMat( SVector(p1, p2, p3, p4) )
    m3_pad = asMatOnePad( SVector(p1, p2, p3) )
    m4_pad = asMatOnePad( SVector(p1, p2, p3, p4) )
    @test all(m3 .== hcat(p1, p2, p3))
    @test all(m4 .== hcat(p1, p2, p3, p4))
    @test all(m3_pad .== hcat(onePad(p1), onePad(p2), onePad(p3)))
    @test all(m4_pad .== hcat(onePad(p1), onePad(p2), onePad(p3), onePad(p4)))
end
