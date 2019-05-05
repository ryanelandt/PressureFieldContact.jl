
# EXAMPLE: Robertson's reaction from pages 3-4 of Solving Ordinary Differential Equations II
#     (Stiff and Differential-Algebraic Problems)

struct RobertsonSystem
    unused_var::Float64  # to make generalization more obvious to readers
    de::Function
    function RobertsonSystem(de::Function)
        return new(1.0, de)
    end
end

function stiff_de!(xx::Vector{T}, x::Vector{T}, s::RobertsonSystem, t::Float64=0.0) where {T}
    c1 = 0.04
    c2 = 1.0e4
    c3 = 3.0e7
    y1 = x[1]
    y2 = x[2]
    y3 = x[3]
    xx[1] = -c1 * y1 + c2 * y2 * y3
    xx[2] =  c1 * y1 - c2 * y2 * y3 - c3 * y2^2
    xx[3] =                           c3 * y2^2
    return nothing
end

my_stiff_struct = RobertsonSystem(stiff_de!)
x0 = [1.0, 0.0, 0.0]  # see Eq. 1.4
initial_h = 1.0e-4

rr = makeRadauIntegrator(my_stiff_struct, x0, 1.0e-16, 2, 1)
update_h!(rr, initial_h)  # force it to start with h = 1.0e-4
rr.rule.s = 2  # force it to start with Radau5
rr.step.h_max = Inf  # don't want timesteps truncated
n_steps = 10

xx = [x0]  # to avoid putting the loop in a function
t = [0.0]  # to avoid putting the loop in a function
for k = 1:n_steps
    x0 = xx[1]
    h, x_final = solveRadau(rr, x0)
    xx[1] = x_final
    t[1] += h
end

println("time: ", t[1])
println("second term: ", xx[1][2])
println("compare to Fig 1.3")

@testset "Robertson's reaction" begin
    @test 3.45e-5 < xx[1][2] < 3.7e-5
    @test t[1] == initial_h * (2^(n_steps) - 1)  # it can at most double h every step
end
