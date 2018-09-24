struct ContactMaterialProperties
    E_effective::Float64
    mu::Float64
    inv_thickness::Float64
    hc_velocity_damping::Float64
    function ContactMaterialProperties(E_effective::Float64, mu::Float64, thickness::Float64, hc_velocity_damping::Float64)
        (1.0e4 < E_effective < 1.0e8) || error("E_effective in unexpected range.")
        (0.0 <= mu < 3.0) || error("mu in unexpected range.")
        (0.001 <= thickness < 1.0) || error("thickness in unexpected range.")
        inv_thickness = 1 / thickness
        (0.3 <= hc_velocity_damping < 5.0) || error("hc_velocity_damping in unexpected range.")
        return new(E_effective, mu, inv_thickness, hc_velocity_damping)
    end
end

struct InertiaMaterialProperties
    thickness::Union{Nothing,Float64}  # if volume mesh is known thickness isn't needed to calculate inertia
    rho::Float64  # rho is always needed to calculate inertia
    function InertiaMaterialProperties(rho, thickness::Union{Nothing,Float64}=nothing)
        if thickness isa Float64
            (0.001 <= thickness < 0.1) || error("thickness in unexpected range.")
        end
        (50.0 <= rho < 2000.0) || error("rho in unexpected range.")
        return new(thickness, rho)
    end
end

struct MaterialProperties
    contact::Union{Nothing,ContactMaterialProperties}
    inertia::Union{Nothing,InertiaMaterialProperties}
    MaterialProperties(contact::ContactMaterialProperties, inertia::InertiaMaterialProperties) = new(contact, inertia)
    MaterialProperties(contact::ContactMaterialProperties) = new(contact, nothing)
    MaterialProperties(inertia::InertiaMaterialProperties) = new(nothing, inertia)
    MaterialProperties() = new(nothing, nothing)
end

function calculateExtrensicCompliance(mat::MaterialProperties)
    (mat.contact == nothing) && return 0.0
    return 1 / (mat.contact.E_effective * mat.contact.inv_thickness)
end
function calcMutualMu(mat_1::MaterialProperties, mat_2::MaterialProperties)
    is_mu_1 = (mat_1.contact != nothing)
    is_mu_2 = (mat_2.contact != nothing)
    if is_mu_1
        return is_mu_2 ? sqrt(mat_1.contact.mu, mat_2.contact.mu) : mat_1.contact.mu
    else
        if is_mu_2
            return mat_2.contact.mu
        else
            error("both mu_1 and mu_2 are nothing")
        end
    end
end