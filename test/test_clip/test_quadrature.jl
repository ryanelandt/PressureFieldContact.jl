
function verifyQuadRule(qr::TriTetQuadRule{N_zeta, N_Point}) where {N_zeta, N_Point}
    (sum(qr.w) ≈ 1.0) || (return false)
    for k = 1:N_Point
        (sum(qr.zeta[k]) ≈ 1.0) || (return false)
    end
    c = zeros(SVector{N_zeta,Float64})
    for k = 1:N_Point
        c += qr.w[k] * qr.zeta[k]
    end
    @test all(c .≈ (1 / N_zeta))
    return true
end

@testset "quadrature" begin
    # Triangle
    for k = 1:5
        @test verifyQuadRule(getTriQuadRule(k))
    end
    qr_2 = getTriQuadRule(2)
    @test all(sort(qr_2.zeta[1]) .== sort(qr_2.zeta[2]) .== sort(qr_2.zeta[3]))

    # Tetrahedron
    for k = 1:5
        @test verifyQuadRule(getTetQuadRule(k))
    end
    @test_throws ErrorException getTriQuadRule((6))
    @test_throws ErrorException getTetQuadRule((6))
end
