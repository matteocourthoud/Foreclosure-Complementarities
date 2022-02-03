#=``
Project: Foreclosure Complementarities
File: Test file
By: Matteo Courthoud
Date: 07/06/2021
=#

include("init.jl")
include("model_lbd.jl")
include("model_privacy.jl")
include("dynamics.jl")
include("solve.jl")
include("postprocess.jl")
disp(x) = Base.print_matrix(stdout, round.(x, digits=4))

# Parameters
alpha = 0.8
sigma = 7

# 1: (109.93 M allocations: 7.184 GiB, 6.15% gc time, 67.46% compilation time)
game = init.model(alpha=alpha, sigma=sigma, policy="baseline")
@time game = solve.solve_game(game)

# 2: (51.37 M allocations: 4.105 GiB, 11.58% gc time, 0.08% compilation time)
game = init.model(alpha=alpha, sigma=sigma, policy="baseline")
@time game = solve.solve_game(game)



# Test
model = "privacy"
a = 0.9
g = 0.0
s = 9
p0 = 1.5

# 3
game = init.model(alpha=a, gamma=g, sigma=s, policy="baseline", modelname=model, p0=p0)
@time game = solve.solve_game(game)

game = init.model(alpha=a, sigma=s, policy="nomergers", modelname=model, p0=p0)
@time game = solve.solve_game(game)

game = init.model(alpha=a, sigma=s, policy="nopredexitpricing", modelname=model, p0=p0)
@time game = solve.solve_game_predatory(game)

postprocess.compute_sumstats(model);
