#=``
Foreclosure Complementarities
By: Matteo Courthoud
Date: 08/06/2021
=#

using BenchmarkTools
include("input/init.jl")
include("input/solve_lbd.jl")
include("input/predatory.jl")
include("input/postprocess.jl")

# Disp function
disp(x) = Base.print_matrix(stdout, Float32.(x))

# Init
precompile(solve_lbd.solve_game, (Main.init.model,))
precompile(predatory.solve_game, (Main.init.model,))
policies = ["baseline", "nolearning", "nomergers", "datasharing", "limitedmergers", "nopredentrypricing", "nopredexitpricing", "nopredentrybundling", "nopredexitbundling"]

# Delete all existing games and issues file
[rm(f) for f in readdir() if f=="issues.txt"]
[rm("output/games/$f") for f in readdir("output/games/") if occursin("json", f)]

# Single parametrization
for sigma = [3, 7]
    for alpha = [0.3, 0.7]
        for policy in policies
            game = init.model(policy=policy, alpha=alpha, gamma=0.0, sigma=sigma);
            if contains(game.policy, "nopred")
                @time game_solved = predatory.solve_game(game);
            else
                @time game_solved = solve_lbd.solve_game(game);
            end
            postprocess.compute_transitions(game_solved)
            postprocess.compute_timelines(game_solved)
        end
    end
end

# Generate comparative statics
for policy in policies
    for sigma = 1:1:10
        for gamma = 0.0:1.0:1.0
            for alpha = 0.0:0.1:1.0
                game = init.model(policy=policy, alpha=alpha, gamma=gamma, sigma=sigma);
                if contains(game.policy, "nopred")
                    @time game_solved = predatory.solve_game(game);
                else
                    @time game_solved = solve_lbd.solve_game(game);
                end
            end
        end
    end
end

# Compute statistics
postprocess.compute_sumstats();
