@testset "exports" begin
    # Ensure that every exported name is actually defined
    for name in names(Radau)
        @test isdefined(Radau, name)
        !isdefined(Radau, name) && println("MISSING: ", name)
    end
end
