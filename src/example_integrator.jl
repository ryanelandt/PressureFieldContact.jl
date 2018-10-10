
function integrate_scenario_radau(mech_scen::MechanismScenario, x::Vector{Float64}; t_final::Float64=1.0, h_max::Float64=0.001, stages::Int64=3)
    n_dof = num_partials(mech_scen)

    (h_max <= 0.0) && error("timestep h_max (set to $h_max) must by positive")
    (t_final < 0.0) && error("t_final (set to $t_final) most be non-negative")
    (length(x) == n_dof) || error("The length of x ($(length(x))) is different from the number of variables in the scenario ($n_dof). Did you forget the state variables for bristle friction?")

    h = h_max
    t_cumulative = 0.0
    rr = makeRadauIntegrator(n_dof, 1.0e-16, mech_scen)
    n_total = ceil(Int64, t_final / h_max)
    data_k_iter = zeros(Int64, n_total)  # positive integer of Newton iterations for Radau coveragence
    data_time = zeros(Float64, n_total)
    data_state = zeros(Float64, n_total, n_dof)
    verify_bristle_ids!(mech_scen, x)  # make sure than any bristle friction variables are used exactly once
    for k = 1:n_total
        is_first = true
        is_converge = false
        while !is_converge
            is_converge, k_iter, res, x_final = solveRadau(rr, x, h, stages, is_first)
            if is_converge
                data_k_iter[k] = k_iter
                x = x_final * 1.0
                h = min(1.2 * h, h_max)
            else
                (h < 1.0e-6) && error("timestep is very small ($h), something is probably wrong")
                h *= 0.5
                is_first = false
            end
        end
        principal_value!(mech_scen, x)
        t_cumulative += h
        data_time[k] = t_cumulative
        data_state[k, :] = x
    end
    return data_time, data_state, data_k_iter
end
