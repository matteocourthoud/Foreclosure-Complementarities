#=``
Project: Foreclosure Complementarities
File: Test file
By: Matteo Courthoud
Date: 07/06/2021
=#

include("init.jl")
include("model_lbd.jl")
include("dynamics.jl")
include("solve.jl")
include("postprocess.jl")
disp(x) = Base.print_matrix(stdout, round.(x, digits=4))

# Parameters
alpha = 0.8
sigma = 7

# 1: 21.313407 seconds (109.93 M allocations: 7.184 GiB, 6.15% gc time, 67.46% compilation time)
game = init.model(alpha=alpha, sigma=sigma, policy="baseline")
@time game = solve.solve_game(game)

# 2: 3.294991 seconds (51.37 M allocations: 4.105 GiB, 11.58% gc time, 0.08% compilation time)
game = init.model(alpha=alpha, sigma=sigma, policy="baseline")
@time game = solve.solve_game(game)
