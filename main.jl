#=``
Foreclosure Complementarities
By: Matteo Courthoud
Date: 07/06/2021
=#


include("input/init.jl")
include("input/lbd.jl")
include("input/predatory.jl")
include("input/postprocess.jl")

# Disp function
function disp(x)
    xx = Float32.(x)
    Base.print_matrix(stdout, xx)
end

# Init
precompile(lbd.solve_game, (Main.init.model,))
precompile(predatory.solve_game, (Main.init.model,))
policies = ["baseline", "nolearning", "nomergers", "datasharing", "limitedmergers", "nopredentrypricing", "nopredexitpricing", "nopredentrybundling", "nopredexitbundling"]

# Solve one game
for policy in policies
    game = init.model(policy=policy);
    if contains(game.policy, "nopred")
        @time game_solved = predatory.solve_game(game);
    else
        @time game_solved = lbd.solve_game(game);
    end
    postprocess.compute_transitions(game_solved)
    postprocess.compute_timelines(game_solved)
end


# Generate comparative statics
for policy in policies
    for sigma = 1:1:9
        for alpha = 0.1:0.1:0.9
                game = init.model(policy=policy, alpha=alpha, sigma=sigma);
            if contains(game.policy, "nopred")
                game_solved = predatory.solve_game(game);
            else
                game_solved = lbd.solve_game(game);
            end
        end
    end
end
# Compute statistics
postprocess.compute_sumstats();
