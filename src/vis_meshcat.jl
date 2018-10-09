
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
