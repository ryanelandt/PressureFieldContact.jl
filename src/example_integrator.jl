
function integrate_scenario_radau(rr::RadauIntegrator{MechanismScenario{NX,NQ,Dual{Nothing,Float64,NC}},NX,NC,NR,NSM},
        mech_scen::MechanismScenario{NX,NQ,Dual{Nothing,Float64,NC}}, x::Vector{Float64};
        t_start::Float64=0.0, t_final::Float64=1.0, max_steps::Int64=1000, is_disp_count::Bool=false) where {NX,NQ,NC,NR,NSM}

    data_state = Matrix{Float64}(x')
    data_time = [t_start]
    return integrate_scenario_radau(rr, mech_scen, data_state, data_time, t_start=t_start, t_final=t_final,
    max_steps=max_steps, is_disp_count=is_disp_count)
end

function integrate_scenario_radau(rr::RadauIntegrator{MechanismScenario{NX,NQ,Dual{Nothing,Float64,NC}},NX,NC,NR,NSM},
        mech_scen::MechanismScenario{NX,NQ,Dual{Nothing,Float64,NC}}, data_state_in::Matrix{Float64},
        data_time_in::Vector{Float64}; t_start::Float64=0.0, t_final::Float64=1.0, max_steps::Int64=1000,
        is_disp_count::Bool=false) where {NX,NQ,NC,NR,NSM}

    n_dof = num_x(mech_scen)
    s1_state, s2_state = size(data_state_in)
    s_time = length(data_time_in)
    (s2_state == n_dof) || error("The length of x ($(length(x))) is different from the number of variables in the scenario ($n_dof). Did you forget the state variables for bristle friction?")
    (t_start < t_final) || error("t_final must be after t_start")
    (s1_state == s_time) || error("data_state and data_time are different lengths")

    i_start = findlast(data_time_in .<= t_start)
    data_state = vcat(data_state_in[1:i_start, :], zeros(Float64, max_steps, n_dof))
    data_time = vcat(data_time_in[1:i_start], zeros(Float64, max_steps))
    t_cumulative = data_time[i_start]
    x = data_state[i_start, :]

    # verify_bristle_ids!(mech_scen, x)  # make sure than any bristle friction variables are used exactly once
    for k = i_start .+ (1:max_steps)
        is_disp_count && println(k)
        h, x = solveRadau(rr, x, t_cumulative)
        principal_value!(mech_scen, x)
        t_cumulative += h
        data_time[k] = t_cumulative
        data_state[k, :] = x
        if t_final < t_cumulative
            return data_time[1:k], data_state[1:k, :]
        end
    end
    return data_time, data_state
end
