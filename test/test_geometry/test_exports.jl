@testset "exports" begin
    # Ensure that every exported name is actually defined
    for name in names(Geometry)
        @test isdefined(Geometry, name)
        !isdefined(Geometry, name) && println("MISSING: ", name)
    end
end
