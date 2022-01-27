"""Network Effects Functions"""

module solve

using JSON3

include("init.jl")
include("model_lbd.jl")
include("dynamics.jl")
include("postprocess.jl")

"""Export game"""
function export_game(game, dist::Float64, iter::Int)
    # Export
    foldername = "data/games/$(game.modelname)/$(game.policy)/"
    mkpath(foldername)
    open("$foldername/$(game.filename).json", "w") do io
        JSON3.write(io, game)
    end

    # If there were issues, save it
    if (dist>100*game.accuracy)
        open("issues.txt","a") do io
           println(io, string(game.filename, " : dist=", dist))
        end
    elseif (max(game.Q[:,end]...)>0.9)
        row = argmax(game.Q[:,end])
        open("issues.txt","a") do io
           println(io, string(game.filename, " : state=", game.S[row, :],
           " prices=", round.(game.P[row, :],digits=1)))
        end
    elseif (iter>400)
        open("issues.txt","a") do io
           println(io, string(game.filename, " : iter=", iter))
        end
    end
end

"""Load game"""
function load_game(game)
    filename = string(game.filename, ".json")
    if string(filename) in readdir("output/games/")
        print("\n\nGame ", game.filename, " already exists!")
        game = init.import_game(filename)
    end
    return game
end

"""Update V"""
function update_V(game, V1, iter)
    # Compute distance
    dist = max(abs.(game.V - V1)...);
    # Update value function
    r = rand() * (iter > 400)
    V = game.V .* r + V1 .* (1-r);
    # Print update
    if iter % 10 == 0
        print("\rIter ", iter, ": ", dist);
    end
    iter += 1;
    # Return
    return V, dist, iter
end

"""Update prices according to the model"""
function update_V_pricing(game, V2, args...)
    if game.modelname == "lbd"
        return model_lbd.update_V_pricing(game, V2)
    elseif game.modelname == "privacy"
        return model_privacy.update_V_pricing(game, V2, args...)
    end
end

"""Solve the dynamic game"""
function solve_game(game)
    # game = load_game(game)

    # Initialize prices to best reply prices with zero value #TODO: compute_W
    game.P = model_lbd.update_P_BR(game, model_lbd.compute_W(game, game.V));

    # Solve game
    print("\n\nSolving ", game.filename, "\n----------------------\n");
    dist = 1;
    iter = 0;
    while (dist>100*game.accuracy) && (iter<1000)
        V5, game.pr_exit = dynamics.update_V_exit(game, game.V);
        V4, game.pr_entry = dynamics.update_V_entry(game, V5);
        V3, game.pr_merger = dynamics.update_V_merger(game, V4);
        V2, game.pr_bundling = dynamics.update_V_bundling(game, V3);
        V1, game.P, game.Q, game.D, game.PI, game.CS = update_V_pricing(game, V2);
        game.V, dist, iter = update_V(game, V1, iter)
    end

    # Compute sumstats
    postprocess.get_sumstats(game);

    # Export game
    export_game(game, dist, iter)

    return game
end

"""Compute index without bundling"""
function compute_idx_nobundling(S::Matrix{Int})

     # Compute new state, removing bundling
     s = copy(S);
     s[:,6] .= 0;

     # Search index of new states
     idx = init.find_states(s, S);
     return idx
end

"""Compute index without learning"""
function compute_idx_nolearning(S::Matrix{Int}, smax)

     # Compute new state, removing bundling
     s = copy(S);
     s12 = s[:,1:2];
     s12[s12 .> 0] .= smax[1];
     s[:,1:2] = s12;

     # Search index of new states
     idx = init.find_states(s, S);
     return idx
end

"""Solve the dynamic game"""
function solve_game_predatory(game)
    # game = load_game(game)

    # Solve game
    print("\n\nSolving ", game.filename, "\n----------------------\n")

    # Make indexes
    idx_nolearning = compute_idx_nolearning(game.S, game.smax)
    idx_nobundling = compute_idx_nobundling(game.S)

    # Initialize value without incentives
    V_noincentives = game.V

    # Initialize prices to best reply prices with zero value #TODO: compute_W
    game.P = model_lbd.update_P_BR(game, model_lbd.compute_W(game, game.V))

    # Iterate until convergence
    dist = 1
    iter = 0
    while (dist>100*game.accuracy) && (iter<1000)

        # Exit: use different probability for non-predatory incentives
        V5, game.pr_exit = dynamics.update_V_exit(game, game.V);
        if contains(game.policy, "nopredexitpricing")
            V_noincentives, _ = dynamics.update_V_exit(game, V_noincentives, game.pr_exit[idx_nolearning]);
        elseif contains(game.policy, "nopredexitbundling")
            V_noincentives, _ = dynamics.update_V_exit(game, V_noincentives, game.pr_exit[idx_nobundling]);
        else
            V_noincentives, _ = dynamics.update_V_exit(game, V_noincentives, game.pr_exit);
        end

        # Entry: use different probability for non-preemptive incentives
        V4, game.pr_entry = dynamics.update_V_entry(game, V5);
        if contains(game.policy, "nopredentrypricing")
            V_noincentives, _ = dynamics.update_V_entry(game, V_noincentives, game.pr_entry[idx_nolearning]);
        elseif contains(game.policy, "nopredentrybundling")
            V_noincentives, _ = dynamics.update_V_entry(game, V_noincentives, game.pr_entry[idx_nobundling]);
        else
            V_noincentives, _ = dynamics.update_V_entry(game, V_noincentives, game.pr_entry);
        end

        # Mergers: use different incentive value for non-predatory bundling
        V3, game.pr_merger = dynamics.update_V_merger(game, V4);
        V_noincentives, _ = dynamics.update_V_merger(game, V_noincentives, V_noincentives, game.pr_merger);

        # Mergers: use different incentive value for non-predatory bundling
        if occursin(r"nopredexitbundling|nopredentrybundling", game.policy)
            V2, game.pr_bundling = dynamics.update_V_bundling(game, V3, V_noincentives);
        else
            V2, game.pr_bundling = dynamics.update_V_bundling(game, V3);
        end
        V_noincentives, _ = dynamics.update_V_bundling(game, V_noincentives, V_noincentives, game.pr_bundling);

        # Pricing: use different incentive value for non-predatory pricing
        if occursin(r"nopredexitpricing|nopredentrypricing", game.policy)
            V1, game.P, game.Q, game.D, game.PI, game.CS = update_V_pricing(game, V3, V_noincentives);
        else
            V1, game.P, game.Q, game.D, game.PI, game.CS = update_V_pricing(game, V3);
        end
        V_noincentives, _, _, _, _, _ = update_V_pricing(game, V_noincentives, V_noincentives, game.P);

        # Compute distance
        game.V, dist, iter = update_V(game, V1, iter)
    end

    # Compute sumstats
    postprocess.get_sumstats(game);

    # Export game
    export_game(game, dist, iter)
    return game
end

"""Replicate alll results"""
function replicate(alphas::Vector, gammas::Vector, sigmas::Vector, policies::Vector{String})

    # Delete all existing games and issues file
    mkpath("data/games/")
    [rm("data/games/$f") for f in readdir("data/games/") if occursin("json", f)]
    [rm(f) for f in readdir() if f=="issues.txt"]


    # Loop over all policies and parameters
    for policy in policies
        for sigma in sigmas
            for gamma in gammas
                for alpha in alphas

                    # Initialize game
                    game = init.model(policy=policy, alpha=alpha, gamma=gamma, sigma=sigma);

                    # Solve game
                    if contains(game.policy, "nopred")
                        game_solved = solve_game_predatory(game);
                    else
                        game_solved = solve_game(game);
                    end

                    # Add details if special game
                    if (sigma in sigmas[[2, end-2]]) && (alpha in alphas[[2, end-2]]) && (gamma == 0.0)
                        postprocess.compute_transitions(game_solved)
                        #postprocess.compute_timelines(game_solved)
                    end

                end
            end
        end
    end

    # Compute statistics
    mkpath("data/sumstats/")
    [rm("data/sumstats/$f") for f in readdir("data/sumstats/") if occursin("csv", f)]
    postprocess.compute_sumstats();
end

end
