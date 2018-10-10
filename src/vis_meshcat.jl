
function MeshCatMechanisms.set_configuration!(mech_scen::MechanismScenario, mvis::MechanismVisualizer, x::Vector{Float64})
    copyto!(mech_scen.float, x)
    set_configuration!(mvis, mech_scen.float.state.q)
    return nothing
end

function set_body_mesh_visual!(mvis::MechanismVisualizer, mech_scen::MechanismScenario, name::String, color)
    return set_body_mesh_visual!(mvis, mech_scen, name, RGBA{Float32}(color...))
end

function set_body_mesh_visual!(mvis::MechanismVisualizer, mech_scen::MechanismScenario, name::String, color::RGBA{Float32})
    body = findbody(mech_scen.float.state.mechanism, name)
    mesh = findMesh(mech_scen.MeshCache, name)
    my_vis_ele = VisualElement(default_frame(body), asHomogenousMesh(mesh), color, Translation(0,0,0))
    setelement!(mvis, my_vis_ele, name)
    return nothing
end

function play_recorded_data(mvis::MechanismVisualizer, mech_scen::MechanismScenario, data_time::Vector{Float64}, data_state::Matrix{Float64};
        dt::Float64=1/60, slowdown::Float64=1.0)

    (size(data_state,1) == length(data_time)) || error("the length of the time vector ($(length(data_time))) needs to be the same as the rows of the state matrix ($(size(data_state,1)))")

    t_last_frame = -dt
    for k = 1:length(data_time)
        if (t_last_frame + dt) < data_time[k]
            t_last_frame += dt
            set_configuration!(mech_scen, mvis, data_state[k, :])
            sleep(dt * slowdown)
        end
    end
    return nothing
end
