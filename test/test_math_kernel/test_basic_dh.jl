
@testset "basic_dh" begin

    zero_3 = zeros(SVector{3,Float64})
    rot_I = one(RotMatrix{3,Float64})
    sm44_I = one(SMatrix{4,4,Float64,16})
    sm33_I = one(SMatrix{3,3,Float64,9})
    sm33_rand = rand(SMatrix{3,3,Float64,9})
    rot_rand = rand(RotMatrix{3,Float64})
    t_rand = rand(SVector{3,Float64})

    # Rotation Only
    @test sm44_I == basic_dh(rot_I).mat
    R, t = dh_R_t(basic_dh(rot_rand))
    @test R == rot_rand
    @test t == zero_3

    # Translation Only
    @test sm44_I == basic_dh(zero_3).mat
    R, t = dh_R_t(basic_dh(t_rand))
    @test R == rot_I
    @test t == t_rand

    # Rotation and Translation
    @test sm44_I == basic_dh(rot_I, zero_3).mat
    R, t = dh_R_t(basic_dh(rot_rand, t_rand))
    @test R == rot_rand
    @test t == t_rand

    # Matrix
    @test sm44_I == basic_dh(sm44_I).mat

    @test R isa SMatrix{3,3,Float64,9}
    @test t isa SVector{3,Float64}

    # SMatrix and Translation
    @test sm44_I == basic_dh(sm33_I, zero_3).mat
    R, t = dh_R_t(basic_dh(sm33_rand, t_rand))
    @test R == sm33_rand
    @test t == t_rand

    # Scale 3
    dh = basic_dh(Diagonal(t_rand))
    @test dh.mat == diagm(0=>[t_rand..., 1.0])

    # Scale constant
    c = 9.0
    dh = basic_dh(c)
    @test dh.mat == diagm(0=>[c, c, c, 1.0])

    # SMatrix only
    @test sm44_I == basic_dh(sm33_I).mat
    R, t = dh_R_t(basic_dh(sm33_rand))
    @test R == sm33_rand
    @test t == zero_3

    # inverse
    dh_rand = basic_dh(rot_rand, t_rand)
    @test inv(dh_rand).mat ≈ inv(dh_rand.mat)

    # multiply
    @test dh_rand.mat * dh_rand.mat ≈ (dh_rand * dh_rand).mat

    # povray_12
    @test all( povray_12(dh_rand) .== SVector(rot_rand..., t_rand...) )

    # one
    @test one(basic_dh{Float64}).mat == zeros(4, 4) + I

    # dh_vector_mul
    p_rand = rand(SVector{3,Float64})
    dh = basic_dh(rot_rand, t_rand)
    @test dh_vector_mul(dh, p_rand) == rot_rand * p_rand + t_rand

end
