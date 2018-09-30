# function addContactRigidCompliant!(mech_scen::MechanismScenario, name_tri::String, name_tet::String)
#     return addContactRigidCompliant!(mech_scen, name_tri, name_tet, nothing)
# end
# function addContactRigidCompliant!(mech_scen::MechanismScenario, name_tri::String, name_tet::String, friction_model::Union{Nothing,BristleFriction})
#     function find_mesh_id(mech_scen::MechanismScenario, name::String)
#         id = MeshID(-9999)
#         for k = mech_scen.mesh_ids
#             if mech_scen.MeshCache[k].name == name
#                 (id == MeshID(-9999)) || error("multiple")
#                 id = k
#             end
#         end
#         return id
#     end
#
#     mesh_id_tri = find_mesh_id(mech_scen, name_tri)
#     mesh_id_tet = find_mesh_id(mech_scen, name_tet)
#     (1 <= mesh_id_tri) || error("invalid tri mesh id $mesh_id_tri")
#     (1 <= mesh_id_tet) || error("invalid tet mesh id $mesh_id_tet")
#     (mesh_id_tri == mesh_id_tet) && error("tri_mesh and tet_mesh id are the same $mesh_id_tri")
#     mesh_cache_tri = mech_scen.MeshCache[mesh_id_tri]
#     mesh_cache_tet = mech_scen.MeshCache[mesh_id_tet]
#
#     (mesh_cache_tet.tet == nothing) && error("compliant mesh named $name_tet has no volume mesh")
#     mat_tet = mesh_cache_tet.tet.contact_prop
#     if mesh_cache_tri.tet == nothing
#         mu = mat_tet.mu
#         frac_epsilon = 1.0
#     else
#         mat_tri = mesh_cache_tri.contact_prop
#         mu = calcMutualMu(mat_tri, mat_tet)
#         hC_tri = calculateExtrensicCompliance(mat_tri)
#         hC_tet = calculateExtrensicCompliance(mat_tet)
#         (hC_tet == 0.0) && error("compliance f tet mesh is rigid because its compliance is zero")
#         frac_epsilon = hC_tet / (hC_tri + hC_tet)
#
#     end
#
#     (0.0 <= mu <= 3.0) || error("mu our of range")
#     frac_linear_weight = 1.0
#     new_contact = ContactInstructions(mesh_id_tri, mesh_id_tet, frac_epsilon, frac_linear_weight, mu, friction_model)
#     push!(mech_scen.ContactInstructions, new_contact)
#     return nothing
# end

# function calc_linear_stiffness(m::MechanismScenario, mesh_id::MeshID; f_disp::Float64=0.0025)
#     gravity = norm(m.float.state.mechanism.gravitational_acceleration.v)
#     mc = mech_scen.MeshCache[mesh_id]
#     char_length = sum(mc.tri.tree.box.e) / 3  # to avoid importing Statistics
#     F = volume(mc) * mc.InertiaProperties.rho
#     delta_x = char_length * f_disp
#     return F / delta_x
# end
#
# function calc_angular_stiffness(m::MechanismScenario, mesh_id::MeshID; rad_disp::Float64=deg2rad(0.25))
#     m = mech_scen
#     mesh_id = MeshID(1)
#     mc = mech_scen.MeshCache[mesh_id]
#     inertia_tensor = bodies(m.float.state.mechanism)[mc.BodyID].inertia
#     avg_inertia = sum(svd(inertia_tensor.moment).S) / 3  # to avoid importing Statistics
#     return avg_inertia / rad_disp
# end
#
# function tune_bristle_friction(m::MechanismScenario, mesh_id::MeshID; tau::Float64=10.0, f_disp::Float64=0.0025, rad_disp::Float64=deg2rad(0.25))
#     K_r = calc_linear_stiffness(m, mesh_id, f_disp=f_disp)
#     K_θ = calc_angular_stiffness(m, mesh_id, rad_disp=rad_disp)
#     return BristleFriction(bristle_id, tau=tau, K_θ=K_θ, K_r=K_r)
# end
