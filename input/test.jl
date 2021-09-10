#=``
Project: Foreclosure Complementarities
File: Test file
By: Matteo Courthoud
Date: 07/06/2021
=#

include("init.jl")
include("solve_lbd.jl")
#include("dynamics.jl")
#include("predatory.jl")
#include("postprocess.jl")
disp(x) = Base.print_matrix(stdout, Float32.(x))

# Solve
game = init.model()
game = solve_lbd.solve_game(game);
game1 = solve_lbd.solve_game(init.model(policy="nolearning"));
game2 = solve_lbd.solve_game(init.model(policy="nomergers"));
#game3 = predatory.solve_game(init.model(policy="nopredexitbundling"));
#game4 = predatory.solve_game(init.model(policy="nopredexitpricing"));
