@testset "exports" begin
    # Ensure that every exported name is actually defined
    for name in names(PressureFieldContact)
        @test isdefined(PressureFieldContact, name)
        !isdefined(PressureFieldContact, name) && println("MISSING: ", name)
    end
end
