
function roll_non_degenerate_tet()
    v1 = randn(SVector{3,Float64})
    v2 = randn(SVector{3,Float64})
    v3 = randn(SVector{3,Float64})
    v4 = randn(SVector{3,Float64})
    if volume(v1, v2, v3, v4) < 0.25
        return roll_non_degenerate_tet()
    else
        return v1, v2, v3, v4
    end
end

function make_4_sided()
    v1, v2, v3, v4 = roll_non_degenerate_tet()
    (0.0 < volume(v1, v2, v3, v4)) || error("bad tet")
    A = asMatOnePad(SVector{4,SVector{3,Float64}}(v1, v2, v3, v4))
    inv_A = inv(A)
    n̂ = normalize(randn(SVector{3,Float64}))
    plane = SMatrix{1,4,Float64,4}(n̂[1], n̂[2], n̂[3], randn())
    c = clip_plane_tet(plane, A)
    if length(c) == 4
        return c
    else
        return make_4_sided()
    end
end

dist_from_plane(plane::SMatrix{1,4,Float64,4}, v::SVector{3,Float64}) = dot(v, unPad(plane)) + plane[4]
project_into_plane(plane::SMatrix{1,4,Float64,4}, v::SVector{3,Float64}) = v - dist_from_plane(plane, v) * unPad(plane)
