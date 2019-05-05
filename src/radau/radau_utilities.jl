
# function makeRadauIntegrator(N::Int64, tol::Float64, de_object::T_object, NC::Int64) where {T_object}
#     return RadauIntegrator{T_object, N, 3, NC}(tol, de_object)
# end
#
# function makeRadauIntegrator(x0::Vector{Float64}, tol::Float64, de_object::T_object, NC::Int64) where {T_object}
#     return makeRadauIntegrator(length(x0), tol, de_object, NC)
# end

function makeRadauIntegrator(de_object::T_object, x_or_N::Union{Int64,Vector{Float64}}, tol::Float64=1.0e-16,
    NR::Int64=2, NC::Int64=6) where {T_object}

    NX = ifelse(x_or_N isa Int64, x_or_N, length(x_or_N))
    NSM = radau_rule_to_stage(NR)
    return RadauIntegrator{NX,NC,NR,NSM}(tol, de_object)
end

function get_X_final(rr::RadauIntegrator{T_object, N, n_stage_max}, table::RadauTable{n_stage}) where {n_stage, n_stage_max, N, T_object}
    return rr.ct.X_stage[n_stage] * 1.0
end

@inline get_exponent(table::RadauTable{n_stage}) where {n_stage} = 1 / (1 + n_stage)
@inline get_h(rr::RadauIntegrator) = rr.step.h

@inline function get_table_from_current_s(rr::RadauIntegrator{T_object, N, n_stage_max}) where {n_stage_max, N, T_object}
    return rr.table[rr.rule.s]
end

function print_exit_flag(rr::RadauIntegrator)
    exit_flag = rr.step.exit_flag
    (exit_flag == -9999) && println("never run")
    (exit_flag == 1) && println("iteration limit exceeded")
    (exit_flag == 2) && println("tolerance failure")
    (exit_flag == 3) && println("newton divergence")
    is_exit_flag_success(rr) && println("step successful")
end

@inline set_exit_flag_success(rr::RadauIntegrator) = (rr.step.exit_flag = 0)
@inline set_exit_flag_fail_iter(rr::RadauIntegrator) = (rr.step.exit_flag = 1)
@inline set_exit_flag_fail_tol(rr::RadauIntegrator) = (rr.step.exit_flag = 2)
@inline set_exit_flag_fail_diverge(rr::RadauIntegrator) = (rr.step.exit_flag = 3)

@inline is_exit_flag_success(rr::RadauIntegrator) = (rr.step.exit_flag == 0)
