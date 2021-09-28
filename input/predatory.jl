"""Network Effects Functions"""

module predatory

include("init.jl")
include("solve_lbd.jl")
include("dynamics.jl")
include("postprocess.jl")

"""Compute index without bundling"""
function compute_idx_nobundling(S::Matrix{Int8})

     # Compute new state, removing bundling
     s = copy(S);
     s[:,5] .= 0;

     # Search index of new states
     idx = init.find_states(s, S);
     return idx
end

"""Compute index without learning"""
function compute_idx_nolearning(S::Matrix{Int8}, smax)

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
function solve_game(game)
    # Check if game exists already
    filename = string(game.filename, ".json")
    if string(filename) in readdir("output/games/")
        #print("\n\nGame ", game.filename, " already exists!")
        #game_solved = init.import_game(filename)
        #return game_solved
    end

    # Solve game
    print("\n\nSolving ", game.filename, "\n----------------------\n")

    # Make indexes
    idx_nolearning = compute_idx_nolearning(game.S, game.smax)
    idx_nobundling = compute_idx_nobundling(game.S)

    # Initialize value without incentives
    V_noincentives = game.V

    # Initialize prices to best reply prices
    P = solve_lbd.update_P_BR(game, solve_lbd.compute_W(game, game.V))
    game.P = solve_lbd.correct_P(game, P)

    # Iterate until convergence
    rate = 1
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
        if occursin(r"nopredexitbundling|nopredentrybundling", game.policy)
            V3, game.pr_merger = dynamics.update_V_merger(game, V4, V_noincentives);
        else
            V3, game.pr_merger = dynamics.update_V_merger(game, V4);
        end
        V_noincentives, _ = dynamics.update_V_merger(game, V_noincentives, V_noincentives, game.pr_merger);

        # Pricing: use different incentive value for non-predatory pricing
        if occursin(r"nopredexitpricing|nopredentrypricing", game.policy)
            V1, game.P, game.Q, game.D, game.PI, game.CS = solve_lbd.update_V(game, V3, V_noincentives);
        else
            V1, game.P, game.Q, game.D, game.PI, game.CS = solve_lbd.update_V(game, V3);
        end
        V_noincentives, _, _, _, _, _ = solve_lbd.update_V(game, V_noincentives, V_noincentives, game.P);

        # Compute distance
        #rate = 0.9*rate + 0.1*max(abs.(game.V - V1)...)/dist
        dist = max(abs.(game.V - V1)...);
        # Update value function
        r = rand() * (iter > 400) # * (rate > 1)
        game.V = game.V .* r + V1 .* (1-r)
        # Print update
        if iter % 10 == 0
            print("\rIter ", iter, ": ", dist)
        end
        iter += 1
    end

    # Print game
    if game.verbose
        sumstats = postprocess.get_sumstats(game);
        print("\n\n", sumstats[:,[7; 10:18]], "\n\n");
        #statestats = postprocess.get_statestats(game);
        #print("\n\n", statestats[:,[3,4,6,7,9,10,11]], "\n\n");
    end

    # Export game
    solve_lbd.export_game(game, dist, iter)
    return game
end

end
