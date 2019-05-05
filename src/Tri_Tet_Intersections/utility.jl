
### Triangle ###
function asMatOnePad(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T}
    A = @SMatrix [
    v1[1] v2[1] v3[1];
    v1[2] v2[2] v3[2];
    v1[3] v2[3] v3[3];
    one(T)  one(T)  one(T)
    ]
    return A
end

function asMat(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T}
    A = @SMatrix [
    v1[1] v2[1] v3[1];
    v1[2] v2[2] v3[2];
    v1[3] v2[3] v3[3]
    ]
    return A
end

for funName in (:asMatOnePad, :asMat)
    @eval begin
        function $funName(sv::SVector{3,SVector{3,T}}) where {T}
            return $funName(sv[1], sv[2], sv[3])
        end
        function $funName(sv::NTuple{3,SVector{3,T}}) where {T}
            return $funName(sv[1], sv[2], sv[3])
        end
    end
end

### Tetrahedron ###
function asMatOnePad(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}, v4::SVector{3,T}) where {T}
    A = @SMatrix [
        v1[1]   v2[1]   v3[1]   v4[1];
        v1[2]   v2[2]   v3[2]   v4[2];
        v1[3]   v2[3]   v3[3]   v4[3];
        one(T)  one(T)  one(T)  one(T)
    ]
    return A
end

function asMat(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}, v4::SVector{3,T}) where {T}
    A = @SMatrix [
        v1[1] v2[1] v3[1] v4[1];
        v1[2] v2[2] v3[2] v4[2];
        v1[3] v2[3] v3[3] v4[3]
    ]
    return A
end

for funName in (:asMatOnePad, :asMat)
    @eval begin
        function $funName(sv::SVector{4,SVector{3,T}}) where {T}
            return $funName(sv[1], sv[2], sv[3], sv[4])
        end
        function $funName(sv::NTuple{4,SVector{3,T}}) where {T}
            return $funName(sv[1], sv[2], sv[3], sv[4])
        end
    end
end
