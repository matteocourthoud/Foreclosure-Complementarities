#=``
Foreclosure Complementarities
By: Matteo Courthoud
Date: 07/06/2021
=#

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

# Delete all existing games
[rm("output/games/$f") for f in readdir("output/games/") if occursin("json", f)]

# Generate comparative statics
for policy in policies
    for sigma = 1:1:9
        for alpha = 0.1:0.1:0.9
                game = init.model(policy=policy, alpha=alpha, sigma=sigma);
            if contains(game.policy, "nopred")
                game_solved = predatory.solve_game(game);
            else
                game_solved = solve_lbd.solve_game(game);
            end
            if (sigma==7.0) & (alpha==0.7)
                postprocess.compute_transitions(game_solved)
                postprocess.compute_timelines(game_solved)
            end
        end

    end
end

# Compute statistics
postprocess.compute_sumstats();
