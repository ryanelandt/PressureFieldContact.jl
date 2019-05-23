@testset "exports" begin
    # Ensure that every exported name is actually defined
    for name in names(MathKernel)
        @test isdefined(MathKernel, name)
        !isdefined(MathKernel, name) && println("MISSING: ", name)
    end
end
