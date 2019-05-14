
@testset "centroid" begin
    p1 = SVector(0.0, 0.0, 0.0)
    p2 = SVector(1.0, 0.0, 0.0)
    p3 = SVector(1.0, 1.0, 0.0)
    p4 = SVector(0.0, 1.0, 0.0)
	n̂ = triangleNormal(p1, p2, p3)

    cent = centroid(poly_eight(4, (p1, p2, p3, p4, p1, p1, p1, p1)), n̂)
    @test cent[1] == 1.0
    @test all(SVector(0.5, 0.5, 0.0) .== cent[2])

    cent = centroid(poly_eight(8, (p1, p2, p3, p4, p1, p1, p1, p1)), n̂)
    @test cent[1] == 1.0
    @test all(SVector(0.5, 0.5, 0.0) .== cent[2])

    cent = centroid(poly_eight(5, (p1, p2, p2, p3, p4, p1, p1, p1)), n̂)
    @test cent[1] == 1.0
    @test all(SVector(0.5, 0.5, 0.0) .== cent[2])

    cent = centroid(poly_eight(3, (p1, p2, p4, p1, p1, p1, p1, p1)), n̂)
    @test cent[1] == 0.5
    @test all(SVector(1/3, 1/3, 0.0) .== cent[2])

    cent = centroid(poly_eight(3, (p1, p2, p2, p1, p1, p1, p1, p1)), n̂)
    @test cent[1] == 0.0
    @test !any(isnan.(cent[2]))
end

@testset "poly_eight" begin
    for i_size = 1:8
        # zero_small_coordinates
        for i_vert = 1:i_size  # each filled vertex
            for i_ind = 1:4
                A = rand(8, 4) .+ 0.5  # all positive numbers
                A[i_vert, i_ind] = (rand() - 0.5) * 3.0e-15  # small number
                tuple_vertex = Tuple(SVector{4,Float64}(A[k, :]) for k = 1:8)
                p8 = poly_eight(i_size, tuple_vertex)
                p8 = zero_small_coordinates(p8)
                for i_vert_check = 1:i_size
                    for i_ind_check = 1:4
                        val_now = p8[i_vert_check][i_ind_check]
                        if (i_ind_check == i_ind) && (i_vert_check == i_vert)
                            @test val_now == 0.0
                        else
                            @test val_now == A[i_vert_check, i_ind_check]
                        end
                    end
                end
            end
        end

        # one_pad_then_mul
        p8 = poly_eight(i_size, Tuple(rand(SVector{3,Float64}) for k = 1:8))
        the_scale = SVector(2.0, 3.0, 4.0)
        m = SArray{Tuple{4,4},Float64,2,16}(diagm(0=>[the_scale..., 0.0]))
        p8_post = one_pad_then_mul(m, p8)

        for i_vert = 1:i_size
            @test all(the_scale .* p8[i_vert] .== p8_post[i_vert][1:3])
        end
    end
end
