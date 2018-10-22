
function asHomogenousMesh(meshCache::MeshCache; color::Union{Nothing, RGBA{Float32}}=nothing)
    vec_Face = Face{3, Int32}.(meshCache.tri.ind)
    vec_Point = Point{3, Float32}.(meshCache.point)
    return HomogenousMesh(vertices=vec_Point, faces=vec_Face, color=color)
end

@RigidBodyDynamics.indextype MeshID
const MeshDict{V} = RigidBodyDynamics.IndexDict{MeshID, Base.OneTo{MeshID}, V}
const MeshCacheDict{V} = RigidBodyDynamics.CacheIndexDict{MeshID, Base.OneTo{MeshID}, V}
Base.@propagate_inbounds Base.getindex(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, key::RigidBody) = d[BodyID(key)]
Base.@propagate_inbounds Base.setindex!(d::RigidBodyDynamics.AbstractIndexDict{MeshID}, value, key::RigidBody) = d[BodyID(key)] = value
