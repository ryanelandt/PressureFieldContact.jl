
@testset "mesh" begin
    eM = eMesh{Tri,Tet}()
    push!(eM.point, SVector{3,Float64}(NaN, NaN, NaN))
    @test_throws ErrorException verify_mesh(eM)
    push!(eM.ϵ, NaN)

    append!(eM, output_eMesh_half_plane())
    mesh_remove_unused_points!(eM)
    @test length(eM.point) == 4

    append!(eM, output_eMesh_half_plane())
    mesh_inplace_rekey!(eM)
    @test length(eM.point) == 4

    @test !isempty(eM)
    @test_throws MethodError empty(eM)
    eM = as_tri_eMesh(eM)
    empty!(eM)
    @test isempty(eM)

    eM = as_tri_eMesh(output_eMesh_box())
    @test area(crop_mesh(eM, SMatrix{1,4,Float64,4}(0.0, 0.0, +1.0, 0.0))) ≈ 12.0
    @test area(crop_mesh(eM, SMatrix{1,4,Float64,4}(0.0, 0.0, -1.0, 0.0))) ≈ 12.0
    @test volume(output_eMesh_box()) ≈ 8.0
end

two = 2.0
@testset "half_plane" begin
    eM_hp_f = output_eMesh_half_plane(two, false)
    eM_hp_t = output_eMesh_half_plane(two, true)
    z_val_f = [k[3] for k = get_point(eM_hp_f)]
    z_val_t = [k[3] for k = get_point(eM_hp_t)]
    @test z_val_f == z_val_t
    @test -two == minimum(z_val_f)
    @test 0.0 == maximum(z_val_f)
    @test 4 == n_point(eM_hp_f)
    for k = 1:3
        @test -0.5 ≈ dot(get_point(eM_hp_f)[k], get_point(eM_hp_f)[mod1(k + 1, 3)])
    end
    @test n_tri(eM_hp_f) == 1
    @test n_tri(eM_hp_t) == 4
    @test n_tet(eM_hp_f) == n_tet(eM_hp_t) == 1
end

@testset "sphere" begin
    for k = 1:4
        eM_sphere = output_eMesh_sphere(two, k)
        n_zero = 0
        @test n_tri(eM_sphere) == n_tet(eM_sphere) == 20 * k^2
        for k = 1:n_point(eM_sphere)
            norm_point = norm(get_point(eM_sphere)[k])
            if norm_point == 0.0
                n_zero += 1
            else
                @test two ≈ norm_point
            end
        end
        @test n_zero == 1
        (k == 1) && @test n_point(eM_sphere) == 1 + 12
        (k == 2) && @test n_point(eM_sphere) == 1 + 12 + 30  # each edge bisected
        (k == 3) && @test n_point(eM_sphere) == 1 + 12 + 2 * 30 + 20  # each edge bisected + point in center
        (k == 4) && @test n_point(eM_sphere) == 1 + 12 + 3 * 30 + 3 * 20  # each edge bisected + two more points in center
    end
end
