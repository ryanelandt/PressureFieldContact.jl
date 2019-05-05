struct RadauVectorCache{n_stage,NX,NC}
    seed::NTuple{NC,ForwardDiff.Partials{NC,Float64}}
    store_float::Vector{Float64}
    store_complex::Vector{ComplexF64}
    neg_J::Matrix{Float64}
    x_dual::Vector{ForwardDiff.Dual{Nothing,Float64,NC}}
    xx_dual::Vector{ForwardDiff.Dual{Nothing,Float64,NC}}
    xx_0::Vector{Float64}
    function RadauVectorCache{n_stage_,NX,NC}() where {n_stage_,NX,NC}
        seed_  = ForwardDiff.construct_seeds(ForwardDiff.Partials{NC,Float64})
        store_float_ = Vector{Float64}(undef, NX)
        store_complex_ = Vector{ComplexF64}(undef, NX)
        neg_J_ = Matrix{Float64}(undef, NX, NX)
        x_dual_ = Vector{ForwardDiff.Dual{Nothing,Float64,NC}}(undef, NX)
        xx_dual_ = Vector{ForwardDiff.Dual{Nothing,Float64,NC}}(undef, NX)
        xx_0_ = Vector{Float64}(undef, NX)
        return new(seed_, store_float_, store_complex_, neg_J_, x_dual_, xx_dual_, xx_0_)
    end
end

struct RadauCacheTuple{n_stage,NX}
    Ew_stage::NTuple{n_stage,Vector{ComplexF64}}
    delta_Z_stage::NTuple{n_stage,Vector{ComplexF64}}
    X_stage::NTuple{n_stage,Vector{Float64}}
    F_X_stage::NTuple{n_stage,Vector{Float64}}
    inv_C_stage::NTuple{n_stage,Matrix{ComplexF64}}
    function RadauCacheTuple{n_stage_,NX}() where {n_stage_,NX}
        Ew_stage_ = Tuple(Vector{ComplexF64}(undef, NX) for _ = 1:n_stage_)
        delta_Z_stage_ = Tuple(Vector{ComplexF64}(undef, NX) for _ = 1:n_stage_)
        X_stage_ = Tuple(Vector{Float64}(undef, NX) for _ = 1:n_stage_)
        F_X_stage_ = Tuple(Vector{Float64}(undef, NX) for _ = 1:n_stage_)
        inv_C_stage_ = Tuple(Matrix{ComplexF64}(undef, NX, NX) for _ = 1:n_stage_)
        return new(Ew_stage_, delta_Z_stage_, X_stage_, F_X_stage_, inv_C_stage_)
    end
end

struct RadauTable{n_stage}
    A::Matrix{Float64}
    c::NTuple{n_stage,Float64}
    λ::Vector{ComplexF64}
    T::Matrix{ComplexF64}
    T⁻¹::Matrix{ComplexF64}
    bi::Matrix{Float64}
    b̂::Vector{Float64}
    b̂_0::Float64
    function RadauTable(n_rule::Int64)
        (0 <= n_rule) || error("radau_rule stage number needs to be positive, but is $n_rule")
        (n_rule <= 6) || error("RadauIIA rule $n_rule with $(radau_rule_to_stage(n_rule)) stages of order $(radau_rule_to_order(n_rule)) not implemented")

        A, b, c, λ, T, T⁻¹, bi, b̂, b̂_0 = load_radau_table_from_file(n_rule)

        n_stages = radau_rule_to_stage(n_rule)
        return new{n_stages}(A, c, λ, T, T⁻¹, bi, b̂, b̂_0)
    end
end

mutable struct RadauStep
    h::Float64
    h⁻¹::Float64
    tol_a::Float64
    tol_r::Float64
    tol_newton::Float64
    is_has_prev_step::Bool
    hᵏ⁻¹::Float64
    x_err_norm::Float64
    x_err_normᵏ⁺¹::Float64
    h_max::Float64
    h_min::Float64
    exit_flag::Int64
    function RadauStep(;h=1.0e-4, tol_a=1.0e-4, tol_r=1.0e-4, tol_newton=1.0e-16)
        return new(h, 1 / h, tol_a, tol_r, tol_newton, false, -9999.0, -9999.0, -9999.0, 0.01, 1.0e-8, -9999)
    end
end

mutable struct RadauRule
    s::Int64
    n_increase_cooldown::Int64
    θ::Float64
    θᵏ⁻¹::Float64
    k_iter::Int64
    k_iter_max::Int64
    Ψ_k::Float64
    max_rule::Int64
    function RadauRule(NR::Int64)
        return new(1, 10, -9999.0, -9999.0, -9999, 15, 9999.0, NR)
    end
end

mutable struct RadauDenseOutput{n_stage}
    is_has_X::Bool
    h::Float64
    s::Int64
    X_stage::NTuple{n_stage,Vector{Float64}}
    function RadauDenseOutput{n_stage}(NX) where {n_stage}
        X_stage = Tuple(Vector{Float64}(undef, NX) for _ = 1:n_stage)
        return new(false, -9999.0, -9999, X_stage)
    end
end

struct RadauIntegrator{T_object,NX,NC,NR,NSM}
    table::NTuple{NR,RadauTable}
    step::RadauStep
    rule::RadauRule
    dense::RadauDenseOutput{NSM}
    cv::RadauVectorCache{NSM,NX,NC}
    ct::RadauCacheTuple{NSM,NX}
    de_object::T_object
    function RadauIntegrator{NX,NC,NR,NSM}(tol::Float64, de_object_::T_object_) where {T_object_,NX,NC,NR,NSM}
        table_ = Tuple([RadauTable(k) for k = 1:NR])
        n_stage_max_ = radau_rule_to_stage(NR)
        cv_ = RadauVectorCache{n_stage_max_,NX,NC}()
        ct_ = RadauCacheTuple{n_stage_max_,NX}()
        radau_step = RadauStep(tol_newton=tol)
        radau_rule = RadauRule(NR)
        dense = RadauDenseOutput{n_stage_max_}(NX)
        return new{T_object_,NX,NC,NR,NSM}(table_, radau_step, radau_rule, dense, cv_, ct_, de_object_)
    end
end
