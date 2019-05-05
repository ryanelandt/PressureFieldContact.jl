
function verify_inplane(plane::SMatrix{1,4,Float64,4}, c::poly_eight{3,Float64}) where {N}
    for k = 1:length(c)
        d_from_plane = dist_from_plane(plane, c[k])
        (1.0e-14 < abs(d_from_plane)) && (return false)
    end
    return true
end

function verify_at_least_2_zero(zeta::SVector{4,Float64})
    n_small = sum(abs.(zeta) .< 1.0e-14)
    return 2 <= n_small
end

@testset "clip_plane_tet" begin
    n_3 = 0
    n_4 = 0
    n_0 = 0
    for k = 1:100
        v1, v2, v3, v4 = roll_non_degenerate_tet()
        @test (0.0 < volume(v1, v2, v3, v4))
        A = asMatOnePad(SVector{4,SVector{3,Float64}}(v1, v2, v3, v4))
        inv_A = inv(A)
        for k = 1:100
            n̂ = normalize(randn(SVector{3,Float64}))
            plane = SMatrix{1,4,Float64,4}(n̂[1], n̂[2], n̂[3], randn())
            c = clip_plane_tet(plane, A)
            d1 = dist_from_plane(plane, v1)
            d2 = dist_from_plane(plane, v2)
            d3 = dist_from_plane(plane, v3)
            d4 = dist_from_plane(plane, v4)
            d_1234 = SVector{4,Float64}(d1, d2, d3, d4)
            if length(c) == 3
                @test (1 == sum(d_1234 .< 0.0)) || (1 == sum(0.0 .< d_1234))
                @test triangleNormal(c[1], c[2], c[3]) ≈ plane[1:3]
                @test verify_at_least_2_zero(inv_A * onePad(c[1]))
                @test verify_at_least_2_zero(inv_A * onePad(c[2]))
                @test verify_at_least_2_zero(inv_A * onePad(c[3]))
                n_3 += 1
            elseif length(c) == 4
                @test 2 == sum(d_1234 .< 0.0)
                @test 2 == sum(0.0 .< d_1234)
                @test triangleNormal(c[1], c[2], c[3]) ≈ plane[1:3]
                @test triangleNormal(c[2], c[3], c[4]) ≈ plane[1:3]
                @test triangleNormal(c[3], c[4], c[1]) ≈ plane[1:3]
                @test triangleNormal(c[4], c[1], c[2]) ≈ plane[1:3]
                @test verify_at_least_2_zero(inv_A * onePad(c[1]))
                @test verify_at_least_2_zero(inv_A * onePad(c[2]))
                @test verify_at_least_2_zero(inv_A * onePad(c[3]))
                @test verify_at_least_2_zero(inv_A * onePad(c[4]))
                n_4 += 1
            else
                @test all(d_1234 .<= 0.0) || all(0.0 .<= d_1234)
                n_0 += 1
            end
            @test verify_inplane(plane, c)
        end
    end
    @test n_0 != 0  # meta-test
    @test n_3 != 0  # meta-test
    @test n_4 != 0  # meta-test
end
