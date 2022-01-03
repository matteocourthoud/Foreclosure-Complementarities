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

alpha = 0.8
sigma = 7

# Compare
game = init.model(alpha=alpha, sigma=sigma, policy="baseline")
@time game = solve.solve_game(game)

game = init.model(alpha=alpha, sigma=sigma, policy="nobundling")
@time game = solve.solve_game(game)

game = init.model(alpha=alpha, sigma=sigma, policy="limitedbundling")
@time game = solve.solve_game(game)

game = init.model(alpha=alpha, sigma=sigma, policy="nopredexitpricing")
@time game = solve.solve_game_predatory(game)

game = init.model(alpha=alpha, sigma=sigma, policy="nopredexitbundling")
@time game = solve.solve_game_predatory(game)

game = init.model(alpha=alpha, sigma=sigma, policy="nopredentrybundling")
@time game = solve.solve_game_predatory(game)

postprocess.compute_sumstats();
