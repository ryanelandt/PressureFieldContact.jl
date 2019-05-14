@testset "exports" begin
    # Ensure that every exported name is actually defined
    for name in names(Clip)
        @test isdefined(Clip, name)
        !isdefined(Clip, name) && println("MISSING: ", name)
    end
end
