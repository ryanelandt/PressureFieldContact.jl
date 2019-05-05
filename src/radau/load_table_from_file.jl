
radau_stage_to_order(n::Int64) = 2 * n - 1
radau_rule_to_stage(n::Int64) = 2 * n - 1
radau_rule_to_order(n::Int64) = radau_rule_to_stage(radau_stage_to_order(n))

get_radau_table_path() = joinpath(first(splitdir(Base.pathof(Radau))), "table")

function load_radau_table_from_file(n_rule::Int64)
    radau_table_path = get_radau_table_path()
    radau_table_path_stage = joinpath(radau_table_path, string(n_rule) * "_rule")

    A      = readdlm(joinpath(radau_table_path_stage, "A.txt"), '\t', Float64, '\n')
    A⁻¹    = readdlm(joinpath(radau_table_path_stage, "inv_A.txt"), '\t', Float64, '\n')
    b      = readdlm(joinpath(radau_table_path_stage, "b.txt"), '\t', Float64, '\n')
    c      = readdlm(joinpath(radau_table_path_stage, "c.txt"), '\t', Float64, '\n')
    λ      = readdlm(joinpath(radau_table_path_stage, "lambda.txt"), '\t', Complex{Float64}, '\n')
    T      = readdlm(joinpath(radau_table_path_stage, "T.txt"), '\t', Complex{Float64}, '\n')
    T⁻¹    = readdlm(joinpath(radau_table_path_stage, "inv_T.txt"), '\t', Complex{Float64}, '\n')
    bi     = readdlm(joinpath(radau_table_path_stage, "bi.txt"), '\t', Float64, '\n')
    b_hat  = readdlm(joinpath(radau_table_path_stage, "b_hat.txt"), '\t', Float64, '\n')

    b = vec(b)
    c = Tuple(vec(c))
    λ = vec(λ)
    b_hat = vec(b_hat)

    b̂ = b_hat[2:end]
    b̂_0 = b_hat[1]

    return A, b, c, λ, T, T⁻¹, bi, b̂, b̂_0
end
