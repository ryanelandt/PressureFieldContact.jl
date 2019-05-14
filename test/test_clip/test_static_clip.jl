
function minimum_area(p::poly_eight{3,Float64}, n̂::SVector{3,Float64}, r::SVector{3,Float64})
    isempty(p) && (return -Inf)
    min_area = Inf
    length_clip = length(p)
    for k = 1:length_clip
        area_k = triangle_area((p[k], p[mod1(k + 1, length_clip)], r), n̂)
        min_area = min(area_k, min_area)
    end
    return min_area
end

@testset "static clip" begin
    tol = 1.0e-13
    n_empty = 0
    n_hits = zeros(MVector{8,Int64})
    for k_tet = 1:500000
        if n_hits[8] <= 2
            r_orig = make_4_sided()

            n̂ = triangleNormal(r_orig[1], r_orig[2], r_orig[3])
            n̂_check = triangleNormal(r_orig[2], r_orig[3], r_orig[4])
            (n̂ ≈ n̂_check) || error("quad points do not lie on a plane")
            plane = SMatrix{1,4,Float64,4}(n̂[1], n̂[2], n̂[3], -dot(n̂, r_orig[1]))

            tet = roll_non_degenerate_tet()
            A = asMatOnePad(tet[1], tet[2], tet[3], tet[4])
            inv_A = inv(A)
            ζ_orig = one_pad_then_mul(inv_A, r_orig)
            ζ_clip = clip_in_tet_coordinates(ζ_orig)
            r_clip = mul_then_un_pad(A, ζ_clip)

            for k = 1:length(r_clip)  # confirm that all vertices on plane
                # TODO: adjust the test to use the origional normal of the clipped shape so this tolerance
                # can be tightened
                @test abs(dist_from_plane(plane, r_clip[k])) < (200 * tol)
            end

            for k = 1:30  # test 30 random points
                r_in_plane = project_into_plane(plane, randn(SVector{3,Float64}))
                min_ζ = minimum(inv_A * onePad(r_in_plane))
                min_area_orig = minimum_area(r_orig, n̂, r_in_plane)  # is_inside origional
                min_area_clip = minimum_area(r_clip, n̂, r_in_plane)  # is_inside clip
                if tol < min_area_clip  # inside clipping
                    point_inside_ζ = -tol < min_ζ
                    point_inside_orig = -tol < min_area_orig
                    @test point_inside_ζ && point_inside_orig
                else
                    point_outside_ζ = min_ζ < tol
                    point_outside_orig = min_area_orig < tol
                    @test point_outside_ζ || point_outside_orig
                end
            end

            if isempty(r_clip)
                n_empty += 1
            else
                n_hits[length(r_clip)] += 1
            end
        end
    end
    @test 2 <= n_hits[8]  # sanity check
    @test 1000 < n_empty  # sanity check
end
