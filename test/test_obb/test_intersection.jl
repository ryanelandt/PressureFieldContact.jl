function outputFaceByNumber(k::Int64)
  (k == 1)  && return SVector{3, Float64}(-1,  0,  0)
  (k == 2)  && return SVector{3, Float64}(+1,  0,  0)
  (k == 3)  && return SVector{3, Float64}( 0, -1,  0)
  (k == 4)  && return SVector{3, Float64}( 0, +1,  0)
  (k == 5)  && return SVector{3, Float64}( 0,  0, -1)
  (k == 6)  && return SVector{3, Float64}( 0,  0, +1)
  error("a cube has 6 faces")
end

function outputCornerByNumber(k::Int64)
  (k == 1) && return SVector{3, Float64}(-1, -1, -1)
  (k == 2) && return SVector{3, Float64}(+1, -1, -1)
  (k == 3) && return SVector{3, Float64}(-1, +1, -1)
  (k == 4) && return SVector{3, Float64}(+1, +1, -1)
  (k == 5) && return SVector{3, Float64}(-1, -1, +1)
  (k == 6) && return SVector{3, Float64}(+1, -1, +1)
  (k == 7) && return SVector{3, Float64}(-1, +1, +1)
  (k == 8) && return SVector{3, Float64}(+1, +1, +1)
  error("a cube has 8 corners")
end

function outputEdgeByNumber(k::Int64)
  (k == 1)  && return SVector{3, Float64}( 0, -1, -1)
  (k == 2)  && return SVector{3, Float64}( 0, +1, -1)
  (k == 3)  && return SVector{3, Float64}( 0, -1, +1)
  (k == 4)  && return SVector{3, Float64}( 0, +1, +1)
  (k == 5)  && return SVector{3, Float64}(-1,  0, -1)
  (k == 6)  && return SVector{3, Float64}(+1,  0, -1)
  (k == 7)  && return SVector{3, Float64}(-1,  0, +1)
  (k == 8)  && return SVector{3, Float64}(+1,  0, +1)
  (k == 9)  && return SVector{3, Float64}(-1, -1,  0)
  (k == 10) && return SVector{3, Float64}(+1, -1,  0)
  (k == 11) && return SVector{3, Float64}(-1, +1,  0)
  (k == 12) && return SVector{3, Float64}(+1, +1,  0)
  error("a cube has 12 edges")
end

function face_corner_test(tt::TT_Cache, dir_1::SVector{3,Float64}, dir_2::SVector{3,Float64})
  aabb_1 = OBB(zero_three, box_1_extent, I3)
  aabb_2 = OBB(zero_three, box_2_extent, I3)
  vec_1 = dir_1 .* aabb_1.e
  vec_2 = dir_2 .* aabb_2.e
  R_total = RotMatrix(rotation_between(dir_2, dir_1))
  total_sep_dist = vec_1 + R_total * vec_2
  update_TT_Cache!(tt, total_sep_dist * (1 - tolerance), R_total)
  is_sect_pos = BB_BB_intersect(tt, aabb_1, aabb_2)
  update_TT_Cache!(tt, total_sep_dist * (1 + tolerance), R_total)
  is_sect_neg = BB_BB_intersect(tt, aabb_1, aabb_2)
  return is_sect_pos, is_sect_neg
end

function edge_edge_test(tt::TT_Cache, dir_1::SVector{3,Float64}, theta_axis::Float64, R_box_2::Rotation{3, Float64})
  aabb_1 = OBB(zero_three, one_three, I3)
  aabb_2 = OBB(zero_three, one_three, I3)
  vec_1 = dir_1 .* aabb_1.e
  R_total = RotMatrix(Quat(AngleAxis(theta_axis, vec_1[1], vec_1[2], vec_1[3])) * R_box_2)
  total_sep_dist = vec_1 * 2
  update_TT_Cache!(tt, total_sep_dist * (1 - tolerance), R_total)
  is_sect_pos = BB_BB_intersect(tt, aabb_1, aabb_2)
  update_TT_Cache!(tt, total_sep_dist * (1 + tolerance), R_total)
  is_sect_neg = BB_BB_intersect(tt, aabb_1, aabb_2)
  return is_sect_pos, is_sect_neg
end

tt = TT_Cache()
box_1_extent = SVector(1.0, 2.0, 3.0)
box_2_extent = SVector(2.1, 2.2, 2.3)
one_three    = SVector(1.0, 1.0, 1.0)
zero_three   = SVector(0.0, 0.0, 0.0)
n_face       = 6
n_corner     = 8
n_edge       = 12
tolerance    = 1.0e-6
I3           = one(SMatrix{3,3,Float64,9})

@testset "OBB OBB Face Face Test" begin
  for k_face = 1:n_face
    for k_corner = 1:n_corner
      dir_face = outputFaceByNumber(k_face)
      dir_corner = outputCornerByNumber(k_corner)
      should_be_positive, should_be_negative = face_corner_test(tt, dir_face, dir_corner)
      @test should_be_positive && !should_be_negative
      should_be_positive, should_be_negative = face_corner_test(tt, dir_corner, dir_face)
      @test should_be_positive && !should_be_negative
    end
  end
end

@testset "OBB OBB Edge Edge Test" begin
  for k_edge_1 = 1:n_face
    for theta_axis = rand(15)*2*pi
      for theta_extra = 0:(pi/2):(2*pi)
        dir_edge = outputEdgeByNumber(k_edge_1)
        should_be_positive, should_be_negative = edge_edge_test(tt, dir_edge, theta_axis, RotX(theta_extra))
        @test should_be_positive && !should_be_negative
        should_be_positive, should_be_negative = edge_edge_test(tt, dir_edge, theta_axis, RotY(theta_extra))
        @test should_be_positive && !should_be_negative
        should_be_positive, should_be_negative = edge_edge_test(tt, dir_edge, theta_axis, RotZ(theta_extra))
        @test should_be_positive && !should_be_negative
      end
    end
  end
end
