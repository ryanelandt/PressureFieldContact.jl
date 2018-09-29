@testset "exports" begin
    # Ensure that every exported name is actually defined
    for name in names(SoftContact)
        @test isdefined(SoftContact, name)
        !isdefined(SoftContact, name) && println("MISSING: ", name)
    end
end
