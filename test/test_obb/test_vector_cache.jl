@testset "VectorCacheInt64" begin
    vc = VectorCache{Int64}()
    @test length(vc) == -9999
    empty!(vc)
    @test isempty(vc)
    @test length(vc) == 0
    @test returnNext(vc) isa Int64

    empty!(vc)
    ind_max_1 = vc.ind_max
    for k = 1:ind_max_1
        addCacheItem!(vc, k)
        @test vc[k] == k
    end
    @test !isempty(vc)

    expand!(vc)
    ind_max_2 = vc.ind_max
    @test length(vc) == vc.ind_fill
    @test ind_max_1 < ind_max_2

    vc.ind_fill = vc.ind_max
    i_next = returnNext(vc)
    ind_max_3 = vc.ind_max
    @test ind_max_2 < ind_max_3
end
