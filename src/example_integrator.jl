
function integrate_scenario_radau(rr::RadauIntegrator{MechanismScenario{NQ,Dual{Nothing,Float64,NC}},NX,NC,NR,NSM};
        t_start::Float64=0.0, t_final::Float64=1.0, max_steps::Int64=1000, is_disp_count::Bool=false) where {NX,NQ,NC,NR,NSM}

    x = get_state(rr.de_object)
    data_state = zeros(max_steps + 1, NX)
    data_time = zeros(max_steps + 1)
    data_state[1, :] = x
    return integrate_scenario_radau(rr, data_state, data_time, t_final=t_final, is_disp_count=is_disp_count)
end

function integrate_scenario_radau(rr::RadauIntegrator{MechanismScenario{NQ,Dual{Nothing,Float64,NC}},NX,NC,NR,NSM},
        data_state::Matrix{Float64}, data_time::Vector{Float64}; t_final::Float64=1.0,
        is_disp_count::Bool=false) where {NX,NQ,NC,NR,NSM}

    mech_scen = rr.de_object
    s1_state, s2_state = size(data_state)
    s_time = length(data_time)
    (s2_state == NX) || error("The length of x ($(length(x))) is different from the number of variables in the scenario ($n_dof). Did you forget the state variables for bristle friction?")
    (s1_state == s_time) || error("data_state and data_time are different lengths")
    t_cumulative = data_time[1]
    x = data_state[1, :]
    for k = 2:s1_state
        is_disp_count && println(k)
        if mech_scen.discrete_controller != nothing
            while mech_scen.discrete_controller.t_last <= t_cumulative
                mech_scen.discrete_controller.t_last += mech_scen.discrete_controller.dt
                mech_scen.discrete_controller.control(x, mech_scen, t_cumulative)
            end
        end
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
