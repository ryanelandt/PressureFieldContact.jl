
@testset "mesh" begin
    eM = eMesh{Tri,Tet}()
    push!(eM.point, SVector{3,Float64}(NaN, NaN, NaN))
    @test_throws ErrorException verify_mesh(eM)
    push!(eM.ϵ, NaN)

    append!(eM, eMesh_half_plane())
    mesh_remove_unused_points!(eM)
    @test length(eM.point) == 4

    append!(eM, eMesh_half_plane())
    mesh_inplace_rekey!(eM)
    @test length(eM.point) == 4

    @test !isempty(eM)
    @test_throws MethodError empty(eM)
    eM = as_tri_eMesh(eM)
    empty!(eM)
    @test isempty(eM)

    eM = as_tri_eMesh(eMesh_box())
    @test area(crop_mesh(eM, SMatrix{1,4,Float64,4}(0.0, 0.0, +1.0, 0.0))) ≈ 12.0
    @test area(crop_mesh(eM, SMatrix{1,4,Float64,4}(0.0, 0.0, -1.0, 0.0))) ≈ 12.0
    @test volume(eMesh_box()) ≈ 8.0
end

@testset "delete_triangles" begin
	random_points = [randn(SVector{3,Float64}) for _ = 1:3]
	i3_tri = SVector{3,Int64}(1,2,3)
	i3_tri_oppose = SVector{3,Int64}(i3_tri[2], i3_tri[3], i3_tri[1])

	eM_delete = eMesh(random_points, [i3_tri, i3_tri_oppose])
	@test 2 == delete_triangles!(eM_delete)
	@test isempty(eM_delete)

	eM_delete = eMesh(random_points, [i3_tri, SVector{3,Int64}(2,3,4)])
	@test 0 == delete_triangles!(eM_delete)
	@test !isempty(eM_delete)

	eM_delete = eMesh(random_points, [i3_tri, i3_tri])
	@test_throws ErrorException delete_triangles!(eM_delete)

	eM_delete = eMesh(random_points, [i3_tri, i3_tri_oppose, i3_tri])
	@test_throws ErrorException delete_triangles!(eM_delete)
end

two = 2.0
@testset "half_plane" begin
    eM_hp_f = eMesh_half_plane(two, false)
    eM_hp_t = eMesh_half_plane(two, true)
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
        eM_sphere = eMesh_sphere(two, k)
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

@testset "polygon circle" begin
	for n = 3:64
		eM_circle = PressureFieldContact.Geometry.create_2D_circle(1.0, n=n)
		b = sin(pi / n)  # distance from center to bisector
		a = sqrt(1.0^2 - b^2)  # distance from top of triangle to point on circle
		my_area = n * a * b
		@test my_area ≈ area(eM_circle)
	end
end

@testset "prism" begin
	thick = 0.1
	p1 = SVector(0.0, 0.0, 0.0)
	p2 = SVector(1.0, 0.0, 0.0)
	p3 = SVector(0.0, 1.0, 0.0)
	eM_tri = eMesh([p1, p2, p3], [SVector(1, 2, 3)])
	eM_prism = extrude_mesh(eM_tri, thick)
	vol_ = volume(eM_prism)
	area_ = area(eM_prism)
	@test vol_ ≈ thick * 0.5
	@test area_ ≈ 1.0 + thick * (2.0 + sqrt(2))
end

@testset "cylinder" begin
	for n = 3:64
		rad = 1.0
		height = 1000.0
		eM_cylinder = eMesh_cylinder(rad, height, n=n)
		eM_circle = PressureFieldContact.Geometry.create_2D_circle(rad, n=n)
		area_circle = area(eM_circle)
		@test volume(eM_cylinder) ≈ height * area_circle
		@test area(eM_cylinder) ≈ 2 * area_circle + height * 2 * n * sin(pi / n)
	end
end
