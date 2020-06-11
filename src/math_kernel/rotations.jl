quatErr(q_err::UnitQuaternion{TD}) where {TD} = (0 <= q_err.w) ? SVector( q_err.x,  q_err.y,  q_err.z) : SVector(-q_err.x, -q_err.y, -q_err.z)

function quatErr(q1::UnitQuaternion{TD}, qRef::UnitQuaternion{TD}) where {TD}
  q_err = q1 * inv(qRef)
  return quatErr(q_err)
end

cheapRV(q::UnitQuaternion) = 2 * quatErr(q)
cheapRV(spq::MRP) = cheapRV(UnitQuaternion(spq))

components(q::UnitQuaternion{T}) where {T} = SVector{4,T}(q.w, q.x, q.y, q.z)
components(spq::MRP{T}) where {T} = SVector{3,T}(spq.x, spq.y, spq.z)
