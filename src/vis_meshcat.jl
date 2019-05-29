
function MeshCatMechanisms.MechanismVisualizer(mech_scen::MechanismScenario, vis)
	return MechanismVisualizer(mech_scen.float.state.mechanism, vis)
end

function MeshCatMechanisms.set_configuration!(mech_scen::MechanismScenario, mvis::MechanismVisualizer)
	set_configuration!(mvis, mech_scen.float.state.q)
end

function set_body_mesh_visual!(mvis::MechanismVisualizer, mech_scen::MechanismScenario, m_id::MeshID, color)
    mesh = mech_scen.MeshCache[m_id]
    body = bodies(mech_scen.float.state.mechanism)[mesh.BodyID]
    return set_body_mesh_visual!(mvis, mech_scen, body, m_id, color)
end

function set_body_mesh_visual!(mvis::MechanismVisualizer, mech_scen::MechanismScenario, body::RigidBody{Float64},
        m_id::MeshID, color)

    color = as_rgba(color)
    mc = mech_scen.MeshCache[m_id]
	mesh = mc.mesh
	mesh_32 = HomogenousMesh_32(mesh)
    my_vis_ele = VisualElement(default_frame(body), mesh_32, color, Translation(0,0,0))
    setelement!(mvis, my_vis_ele, mc.name)
    return nothing
end

function set_mesh_visual!(mvis::MechanismVisualizer, mech_scen::MechanismScenario, m_id::MeshID, color)
    mc = mech_scen.MeshCache[m_id]
    body = root_body(mech_scen.float.state.mechanism)
    color = as_rgba(color)
	mesh = mc.mesh
	mesh_32 = HomogenousMesh_32(mesh)
    my_vis_ele = VisualElement(default_frame(body), mesh_32, color, Translation(0,0,0))
    setelement!(mvis, my_vis_ele, mc.name)
    return nothing
end

as_rgba(color::RGBA{Float32}) = color
as_rgba(color) = RGBA{Float32}(color...)

get_vertices_32(e_mesh::eMesh{T1,T2})  where {T1,T2} = [Point{3,Float32}(k) for k = e_mesh.point]
get_faces_32(   e_mesh::eMesh{Tri,T2}) where {T2}    = [Face{3,Int32}(k)    for k = e_mesh.tri]

function HomogenousMesh_32(eM::eMesh)
	eM = as_tri_eMesh(eM, false)
	vertices = get_vertices_32(eM)
    faces = get_faces_32(eM)
    return HomogenousMesh(vertices=vertices, faces=faces)
end

# function HomogenousMesh_32(e_mesh::eMesh{Tri,T2}; color=RGBA{Float32}(0.5, 0.5, 0.5, 1.0)) where {T2}
#     vertices = get_vertices_32(e_mesh)
#     faces = get_faces_32(e_mesh)
#     return HomogenousMesh(vertices=vertices, faces=faces, color=color)
# end
#
# function HomogenousMesh_32(eM::eMesh{Nothing,Tet})
# 	i3 = Vector{SVector{3,Int32}}()
# 	for k = 1:n_tet(eM)
# 		i_tet = eM.tet[k]
# 		ϵ_tet = eM.ϵ[i_tet]
# 		i_tet_new = sort_so_big_ϵ_last(ϵ_tet, i_tet)
# 		push!(i3, i_tet_new[1:3])
# 	end
# 	vertices=get_vertices_32(eM)
# 	faces = [Face{3,Int32}(k) for k = i3]
# 	return HomogenousMesh(vertices=vertices, faces=faces)
# end

function play_recorded_data(mvis::MechanismVisualizer, mech_scen::MechanismScenario, data_time::Vector{Float64},
        data_state::Matrix{Float64}; dt::Float64=1/60, slowdown::Float64=1.0, t0::Float64=-Inf, t1::Float64=Inf)

    ### ERROR checking
    (size(data_state,1) == length(data_time)) || error("the length of the time vector ($(length(data_time))) needs to be the same as the rows of the state matrix ($(size(data_state,1)))")
    issorted(data_time) || error("data_time (vector of times) is not sorted")
    @assert(t0 < t1)
    i0 = findfirst(t0 .<= data_time)
    i1 = findlast(data_time .<= t1)
    delta_t_sim = data_time[i1] - data_time[i0]
    (60.0 < (slowdown * delta_t_sim) ) && error("total video time is greater than 1 minute you probably didn't mean to do this.")

    t_last_frame = data_time[i0] - dt
    for k = i0:i1
        while (t_last_frame + dt) < data_time[k]  # wait until enough real time has passed
            t_last_frame += dt
            sleep(dt * slowdown)
			copyto!(mech_scen.float, data_state[k, :])
			set_configuration!(mech_scen, mvis)
        end
    end
end
