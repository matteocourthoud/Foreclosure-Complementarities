#=``
Project: Foreclosure Complementarities
File: Test file
By: Matteo Courthoud
Date: 07/06/2021
=#

include("input/init.jl")
include("input/lbd.jl")
include("input/dynamics.jl")
include("input/predatory.jl")
include("input/postprocess.jl")
disp(x) = Base.print_matrix(stdout, Float32.(x))

# Params
alpha = 0.7
gamma = 1.0
sigma = 7.0
v0 = 1

# Solve
game = lbd.solve_game(init.model(alpha=alpha, gamma=gamma, sigma=sigma, v0=v0));
game1 = lbd.solve_game(init.model(policy="nolearning", alpha=alpha, gamma=gamma, sigma=sigma, v0=v0));
game1 = lbd.solve_game(init.model(policy="nomergers", alpha=alpha, gamma=gamma, sigma=sigma, v0=v0));
game2 = predatory.solve_game(init.model(policy="nopredexitbundling", alpha=alpha, gamma=gamma, sigma=sigma, v0=v0));
game3 = predatory.solve_game(init.model(policy="nopredexitpricing", alpha=alpha, gamma=gamma, sigma=sigma, v0=v0));
