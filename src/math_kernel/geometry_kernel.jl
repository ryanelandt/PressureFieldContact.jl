
### TRIANGLE
area(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T} = norm(vector_area(v1, v2, v3))
centroid(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T} = (v1 + v2 + v3) * Float64(1/3)  # 4 times faster than dividing by 3
vector_area(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T} = cross(v2 - v1, v3 - v2) * 0.5
function triangle_area(v123::NTuple{3,SVector{3,T}}, n̂::SVector{3,T}) where {T}
	vec_area = vector_area(v123[1], v123[2], v123[3])
	return dot(n̂, vec_area)
end
triangleNormal(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T} = normalize(vector_area(v1, v2, v3))

for funName in (:area, :centroid, :triangleNormal, :triangleCross)
    @eval begin
        function $funName(sv::SVector{3,SVector{3,T}}) where {T}
            return $funName(sv[1], sv[2], sv[3])
        end
    end
end

### TETRAHEDRON
function centroid(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}, v4::SVector{3,T}) where {T}
    return (v1 + v2 + v3 + v4) * 0.25
end

function volume(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}, v4::SVector{3,T}) where {T}
    # NOTE: This is an algebraic refactorization of a symbolic answer
    a1, a2, a3 = v1[1], v1[2], v1[3]
    b1, b2, b3 = v2[1], v2[2], v2[3]
    c1, c2, c3 = v3[1], v3[2], v3[3]
    d1, d2, d3 = v4[1], v4[2], v4[3]
    V =       (b1 - a1) * (c2 * d3 - c3 * d2)
    V = muladd(b2 - a2,    c3 * d1 - c1 * d3, V)
    V = muladd(b3 - a3,    c1 * d2 - c2 * d1, V)
    V = muladd(c1 - d1,    a3 * b2 - a2 * b3, V)
    V = muladd(c2 - d2,    a1 * b3 - a3 * b1, V)
    V = muladd(c3 - d3,    a2 * b1 - a1 * b2, V)
    return V * Float64(1/6)
end

for funName in (:centroid, :volume)
    @eval begin
        function $funName(sv::SVector{4,SVector{3,T}}) where {T}
            return $funName(sv[1], sv[2], sv[3], sv[4])
        end
    end
end
