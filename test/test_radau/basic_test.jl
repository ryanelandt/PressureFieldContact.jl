struct SimpleStiffSystem
    λ::Float64
    de::Function
    function SimpleStiffSystem(de::Function)
        return new(1.0, de)
    end
end

function stiff_de!(xx::Vector{T}, x::Vector{T}, s::SimpleStiffSystem, t::Float64=0.0) where {T}
    xx .= -s.λ * x
    return nothing
end

my_stiff_struct = SimpleStiffSystem(stiff_de!)
x0 = zeros(Float64, 4) .+ 1.0
n_rule_max = 3

@testset "basic_test" begin
    for NC = 1:10  # check different chunk sizes
        rr = makeRadauIntegrator(my_stiff_struct, x0, 1.0e-16, n_rule_max, NC)
        for k_rule = 1:n_rule_max  # chick different rules
            rr.rule.s = 3
            t_final = 0.2
            update_h!(rr, t_final)
            h, x_final = solveRadau(rr, x0)
            if k_rule != 1
                x_ana = exp(-t_final)
                @test x_ana ≈ x_final[1]
            end
        end
    end
end
