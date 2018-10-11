struct ContactProperties
    Ē::Float64
    μ::Float64
    d⁻¹::Float64
    χ::Float64
    function ContactProperties(;Ē::Float64, μ::Float64, d::Float64, χ::Float64)
        (1.0e4 <= Ē <= 1.0e8) || error("E_effective in unexpected range.")
        (0.0 <= μ <= 3.0) || error("mu in unexpected range.")
        (0.001 <= d <= 1.0) || error("thickness in unexpected range.")
        d⁻¹ = 1 / d
        (0.3 <= χ <= 5.0) || error("hc_velocity_damping in unexpected range.")
        return new(Ē, μ, d⁻¹, χ)
    end
end

struct InertiaProperties
    d::Union{Nothing,Float64}  # if volume mesh is known thickness isn't needed to calculate inertia
    rho::Float64  # rho is always needed to calculate inertia
    function InertiaProperties(;rho::Union{Nothing,Float64}=nothing, d::Union{Nothing,Float64}=nothing)
        (rho == nothing) && error("rho is required")
        if d isa Float64
            (0.001 <= d < 0.1) || error("thickness in unexpected range.")
        end
        (50.0 <= rho < 2000.0) || error("rho in unexpected range.")
        return new(d, rho)
    end
end

function calculateExtrensicCompliance(mat::ContactProperties)
    return 1 / (mat.Ē * mat.d⁻¹)
end

function calcMutualMu(mat_1::ContactProperties, mat_2::ContactProperties)
    error("does anything actually use this function?")
    if is_mu_1
        return is_mu_2 ? sqrt(mat_1.mu, mat_2.mu) : mat_1.mu
    else
        if is_mu_2
            return mat_2.mu
        else
            error("both mu_1 and mu_2 are nothing")
        end
    end
end
