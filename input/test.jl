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
disp(x) = Base.print_matrix(stdout, round.(x, digits=4))


# Having a complementary market increases predatory pricing and market monopolization
game1 = init.model(alpha=0.4, gamma=0.0, sigma=7.0, mergers=false)
@time game1 = solve_lbd.solve_game(game1);

game2 = init.model(alpha=0.4, gamma=1.0, sigma=7.0, mergers=false)
@time game2 = solve_lbd.solve_game(game2);

# How much of it is predatory? All of it
game1 = init.model(alpha=0.4, gamma=0.0, sigma=7.0, mergers=false, policy="nopredexitpricing")
@time game1 = predatory.solve_game(game1);

game2 = init.model(alpha=0.4, gamma=1.0, sigma=7.0, mergers=false, policy="nopredexitpricing")
@time game2 = predatory.solve_game(game2);

# What if now firms can also merge and bundle? Even worse
game1 = init.model(alpha=0.4, gamma=0.0, sigma=7.0)
@time game1 = solve_lbd.solve_game(game1);

game2 = init.model(alpha=0.4, gamma=1.0, sigma=7.0)
@time game2 = solve_lbd.solve_game(game2);

# How much of it is predatory?
game1 = init.model(alpha=0.4, gamma=0.0, sigma=7.0, policy="nopredexitpricing")
@time game1 = predatory.solve_game(game1);

game2 = init.model(alpha=0.4, gamma=1.0, sigma=7.0, policy="nopredexitpricing")
@time game2 = predatory.solve_game(game2);

# How much of it is predatory?
game1 = init.model(alpha=0.4, gamma=0.0, sigma=7.0, policy="nopredexitbundling")
@time game1 = predatory.solve_game(game1);

game2 = init.model(alpha=0.4, gamma=1.0, sigma=7.0, policy="nopredexitbundling")
@time game2 = predatory.solve_game(game2);


# Test
game1 = init.model(alpha=0.4, gamma=0.0, sigma=7.0)
@time game1 = solve_lbd.solve_game(game1);

game1 = init.model(alpha=0.4, gamma=0.0, sigma=7.0, mergers=false)
@time game1 = solve_lbd.solve_game(game1);

game1 = init.model(alpha=0.1, gamma=0.0, sigma=8.0)
@time game1 = solve_lbd.solve_game(game1);

game1 = init.model(alpha=0.1, gamma=0.0, sigma=8.0, mergers=false)
@time game1 = solve_lbd.solve_game(game1);



include("init.jl")
for sigma = 1:1:10
    for alpha = 0.0
        game = init.model(alpha=alpha, gamma=0.0, sigma=sigma);
        @time game_solved = solve_lbd.solve_game(game);
    end
end


game = init.model(alpha=0.4, gamma=0.0, sigma=8.0, policy="nobundling")
@time game = solve_lbd.solve_game(game);

game = init.model(alpha=0.4, gamma=0.0, sigma=8.0, policy="nomergers")
@time game = solve_lbd.solve_game(game);
