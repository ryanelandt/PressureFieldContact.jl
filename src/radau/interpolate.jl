
function interpolate_output!(x_out::Vector{Float64}, rr::RadauIntegrator{T_object, N, n_stage_max}, h::Float64, is_extrapolate_prev::Bool=true) where {n_stage_max, N<:Int64, T_object}
    rr.dense.is_has_X || error("no dense output to interpolate")
    dense = rr.dense
    theta = h / dense.h + is_extrapolate_prev
    theta_prod = theta
    n_stage = dense.s
    theta_vec = ones(n_stage)
    theta_vec = cumprod(theta_vec)
    table = rr.table[n_stage]
    bi = table.bi
    weights = bi * theta_vec
    LinearAlgebra.BLAS.blascopy!(N, dense.X_stage[1], 1, x_out, 1)
    weight_1 = weights[1]
    LinearAlgebra.BLAS.scal!(N, weight_1, x_out, 1)
    for k = 2:n_stage
        weight_k = weights[k]
        LinearAlgebra.BLAS.axpy!(weight_k, dense.X_stage[k], x_out)
    end
    return nothing
end

function initialize_X_with_interp!(rr::RadauIntegrator{T_object, N, n_stage_max}, table::RadauTable{n_stage}) where {n_stage, n_stage_max, N, T_object}
    h = rr.step.h
    c = table.c
    for k = 1:n_stage
        interpolate_output!(rr.ct.X_stage[k], rr, h * c[k], true)
    end
    return nothing
end

@inline function update_dense_unsuccessful!(rr::RadauIntegrator{T_object, N, n_stage_max}) where {n_stage_max, N, T_object}
    rr.dense.is_has_X = false
    return nothing
end

function update_dense_successful!(rr::RadauIntegrator{T_object, N, n_stage_max}, table::RadauTable{n_stage}) where {n_stage_max, N, T_object, n_stage}
    rr.dense.is_has_X = true
    rr.dense.h = rr.step.h
    rr.dense.s = n_stage
    for i = 1:n_stage
        LinearAlgebra.BLAS.blascopy!(N, rr.ct.X_stage[i], 1, rr.dense.X_stage[i], 1)
    end
    return nothing
end
