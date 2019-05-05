@testset "exports" begin
    # Ensure that every exported name is actually defined
    for name in names(Binary_BB_Trees)
        @test isdefined(Binary_BB_Trees, name)
        !isdefined(Binary_BB_Trees, name) && println("MISSING: ", name)
    end
end
