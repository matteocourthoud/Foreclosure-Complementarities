#=``
Project: Foreclosure Complementarities
File: Test file
By: Matteo Courthoud
Date: 07/06/2021
=#

include("init.jl")
include("solve_lbd.jl")
include("dynamics.jl")
include("predatory.jl")
include("postprocess.jl")
disp(x) = Base.print_matrix(stdout, Float32.(x))



game1 = init.model(alpha=0.4, gamma=0.0, sigma=7.0, mergers=false, policy="nopredexitpricing")
@time game1 = predatory.solve_game(game1);

game2 = init.model(alpha=0.4, gamma=1.0, sigma=7.0, mergers=false, policy="nopredexitpricing")
@time game2 = predatory.solve_game(game2);

# Solve
game1 = init.model(alpha=0.4, gamma=0.0, sigma=7.0)
@time game1 = solve_lbd.solve_game(game1);

game2 = init.model(alpha=0.4, gamma=1.0, sigma=7.0)
@time game2 = solve_lbd.solve_game(game2);

game1 = init.model(alpha=0.4, gamma=0.0, sigma=7.0, policy="nopredexitpricing")
@time game1 = predatory.solve_game(game1);

game2 = init.model(alpha=0.4, gamma=1.0, sigma=7.0, policy="nopredexitpricing")
@time game2 = predatory.solve_game(game2);




@time game1 = solve_lbd.solve_game(init.model(policy="nolearning"));
@time game2 = solve_lbd.solve_game(init.model(policy="nomergers"));
@time game3 = predatory.solve_game(init.model(policy="nopredexitbundling"));
@time game4 = predatory.solve_game(init.model(policy="nopredexitpricing"));



game1 = init.model(alpha=0.4, gamma=0.0, sigma=7.0)
@time game1 = solve_lbd.solve_game(game1);

game1 = init.model(alpha=0.4, gamma=0.5, sigma=7.0)
@time game1 = solve_lbd.solve_game(game1);


# Solve
game1 = init.model(alpha=0.4, gamma=0.0, sigma=2.0)
@time game1 = solve_lbd.solve_game(game1);

game2 = init.model(alpha=0.4, gamma=0.0, sigma=2.0, policy="nomergers")
@time game2 = solve_lbd.solve_game(game2);
