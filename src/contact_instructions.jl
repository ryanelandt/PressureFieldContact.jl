# addContactRigidCompliant!

function addContactRigidCompliant!(mech_scen::MechanismScenario, name_tri::String, name_tet::String)
    return addContactRigidCompliant!(mech_scen, name_tri, name_tet, nothing)
end
function addContactRigidCompliant!(mech_scen::MechanismScenario, name_tri::String, name_tet::String, friction_model::Union{Nothing,BristleFriction})
    function find_mesh_id(mech_scen::MechanismScenario, name::String)
        id = MeshID(-9999)
        for k = mech_scen.mesh_ids
            if mech_scen.MeshCache[k].name == name
                (id == MeshID(-9999)) || error("multiple")
                id = k
            end
        end
        return id
    end

    mesh_id_tri = find_mesh_id(mech_scen, name_tri)
    mesh_id_tet = find_mesh_id(mech_scen, name_tet)
    (1 <= mesh_id_tri) || error("invalid tri mesh id $mesh_id_tri")
    (1 <= mesh_id_tet) || error("invalid tet mesh id $mesh_id_tet")
    (mesh_id_tri == mesh_id_tet) && error("tri_mesh and tet_mesh id are the same $mesh_id_tri")
    mesh_cache_tri = mech_scen.MeshCache[mesh_id_tri]
    mesh_cache_tet = mech_scen.MeshCache[mesh_id_tet]
    mat_tet = mesh_cache_tet.raw.material
    mat_tri = mesh_cache_tri.raw.material
    if mat_tet == nothing
        error("tet mesh must be compliant but there are no MaterialProperties")
    elseif mat_tri == nothing
        mu = mat_tet.contact.mu
        frac_epsilon = 1.0
    else
        mu = calcMutualMu(mat_tri, mat_tet)
        hC_tri = calculateExtrensicCompliance(mat_tri)
        hC_tet = calculateExtrensicCompliance(mat_tet)
        (hC_tet == 0.0) && error("compliance f tet mesh is rigid because its compliance is zero")
        frac_epsilon = hC_tet / (hC_tri + hC_tet)
    end
    (0.0 <= mu <= 3.0) || error("mu our of range")
    frac_linear_weight = 1.0
    new_contact = ContactInstructions(mesh_id_tri, mesh_id_tet, frac_epsilon, frac_linear_weight, mu, friction_model)
    push!(mech_scen.ContactInstructions, new_contact)
    return nothing
end
