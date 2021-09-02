"""Entry, exit, mergers and bundling functions"""

module dynamics

"""Compute probability"""
function compute_p(Delta::Vector{Float64}, c::Vector{Int8})::Tuple{Vector{Float64},Vector{Float64}}
    p = min.(max.( (Delta .- c[1]) ./ (c[2]-c[1]), 0), 1) ./ 4;
    ev = max.(c[1] .+ (min.(Delta, c[2]) .- c[1]) ./ 2, 0);
    return p, ev
end

"""Update value for exit"""
function update_V_exit(game, V::Matrix{Float64}, args...)::Tuple{Matrix{Float64},Matrix{Float64}}

    # Init
    V_exit = Float64.(zeros(size(V)));
    pr_exit = Float64.(zeros(size(game.pr_exit)));

    # Return if no exit
    if game.exit == false
        return V, pr_exit
    end

    # Calculate exit probability and update value for each firm
    for e=1:4

        # Compute exit value
        V1 = V[game.idx_exit[:,e,:]];

        # Compute delta
        Delta = V[:,e] - V1[:,e];

        # Correct Delta: take into account effect of own exit on partner
        rows = (game.ownership[game.S[:,5] .+ 1, e] .> 0)
        Delta_partner = V[rows,game.partner[e]] - V1[rows,game.partner[e]];
        Delta[rows] += Delta_partner;

        # Compute policy if not provided
        if length(args) == 0
            # Compute exit probability
            v = game.value_exit;
            p, _ = compute_p(Delta, v);
            pr_exit[:,e] = (1 .- 4 .* p) ./ 4; # convert pr(non-exit) into pr(exit)

            # Cannot exit if already out or if alone in the market
            no_exit = (game.S[:,e] .== 0) + (game.S[:,game.rival[e]] .== 0);
            pr_exit[no_exit .> 0, e] .= 0;
        else
            pr_exit[:,e] = args[1][:,e]
        end

        # Update value of each other firm
        V1[:,e] = V[:,e];
        V_exit += V1 .* pr_exit[:,e];
    end

    # Update value with residual probability
    V_exit += V .* (1 .- sum(pr_exit, dims=2));

    @assert max(abs.(V_exit[game.S[:,1:4].==0])...) == 0
    return V_exit, pr_exit
end

"""Update value for entry"""
function update_V_entry(game, V::Matrix{Float64}, args...)::Tuple{Matrix{Float64},Matrix{Float64}}

    # Init value
    V_entry = Float64.(zeros(size(V)));
    pr_entry = Float64.(zeros(size(game.pr_entry)));

    # Return if no entry
    if game.entry == false
        return V, pr_entry
    end

    # Calculate entry probability and update value for each firm
    for e=1:4

        # Compute entry value
        V1 = V[game.idx_entry[:,e,:]];

        # Compute entry surplus
        Delta = V1[:,e] - V[:,e];

        # Compute policy if not provided
        c = game.cost_entry;
        if length(args) == 0
            pr_entry[:,e], exp_cost = compute_p(Delta, c);
        else # CHECK but should be correct
            pr_entry[:,e] = args[1][:,e]
            exp_cost = (c[1] .+ (c[2] .* pr_entry[:,e] + c[1] .* (1 .- pr_entry[:,e])) ./2 ) .* pr_entry[:,e]
        end

        # Update entry value
        V1[:,e] -= exp_cost;
        V1[game.S[:,1:4] .== 0] .= 0;
        V_entry += V1 .* pr_entry[:,e];
    end

    # Update value with residual probability
    V_entry += V .* (1 .- sum(pr_entry, dims=2));

    @assert max(abs.(V_entry[game.S[:,1:4].==0])...) == 0
    return V_entry, pr_entry
end

"""Update value for mergers"""
function update_V_merger(game, V::Matrix{Float64}, args...)::Tuple{Matrix{Float64},Matrix{Float64}}

    # Init
    V_merger = Float64.(zeros(size(V)));
    pr_merger = Float64.(zeros(size(game.pr_merger)));

    # Return if no mergers
    if game.mergers == false
        return V, pr_merger
    end

    # Update with incentives
    V_incentives = (length(args)==1) ? args[1] : V

    # Loop over all possible merger pairs
    for p=1:size(game.merger_pairs,1)

        # Compute merger value
        V1 = V[game.idx_merger[:,p,:]];
        V1_incentives = V_incentives[game.idx_merger[:,p,:]];

        # Compute policy if not provided
        c = game.cost_merger;
        Delta = reshape(sum(V1[:,game.merger_pairs[p,:]] - V[:,game.merger_pairs[p,:]], dims=2), game.ms);
        Delta_incentives = reshape(sum(V1_incentives[:,game.merger_pairs[p,:]] - V_incentives[:,game.merger_pairs[p,:]], dims=2), game.ms);
        if length(args) < 2
            pr_merger[:,p], exp_cost = compute_p(Delta_incentives, c);
        else # CHECK but should be correct
            pr_merger[:,p] = args[2][:,p]
            exp_cost = (c[1] .+ (c[2] .* pr_merger[:,p] + c[1] .* (1 .- pr_merger[:,p])) ./2 ) .* pr_merger[:,p]
        end

        # Update value with outcome of Nash bargaining: CHECK
        V1[:,game.merger_pairs[p,:]] = V[:,game.merger_pairs[p,:]] .+ (Delta-exp_cost) ./ 2
        V_merger += V1 .* pr_merger[:,p];
    end

    # Add residual value
    V_merger += V .* (1 .- sum(pr_merger, dims=2));

    @assert max(abs.(V_merger[game.S[:,1:4].==0])...) == 0
    return V_merger, pr_merger
end


"""Update value for bundling"""
function update_V_bundl(game, V::Matrix{Float64}, args...)::Tuple{Matrix{Float64},Matrix{Float64}}

    # Init
    V_bundl = Float64.(zeros(size(V)));
    pr_bundl = Float64.(zeros(size(game.pr_bundl)));

    # Return if no bundling
    if game.bundling == false
        return V, pr_bundl
    end

    # Init present and future value
    V1 = V[game.idx_bundl];
    V_incentives = (length(args)==1) ? args[1] : V
    V1_incentives = V_incentives[game.idx_bundl];

    # Init policy and cost
    exp_cost = Float64.(zeros(game.ms,4));

    # Calculate bundling probability and update value for each firm
    for n=1:4

        # Rows where firms can decide to bundle
        ownership = game.ownership[game.S[:,5] .+ 1, n]
        rows = ((ownership .> 0) .* (game.S[:,6] .== 0)) .> 0;

        # Compute Delta: difference in value with and without bundling
        Delta = Float64.(zeros(game.ms));
        Delta[rows] = V1_incentives[rows,n] + V1_incentives[rows,game.partner[n]] - V_incentives[rows,n] - V_incentives[rows,game.partner[n]];

        # Compute policy if not provided
        c = game.cost_bundl;
        if length(args) < 2
            pr_bundl[:,n], exp_cost[:,n] = compute_p(Delta, c);
        else # CHECK but should be correct
            pr_bundl[:,n] = args[2][:,n]
            exp_cost[:,n] = (c[1] .+ (c[2] .* pr_bundl[:,n] + c[1] .* (1 .- pr_bundl[:,n])) ./2 ) .* pr_bundl[:,n]
        end
    end

    # Update V
    V_bundl = (V1 .- exp_cost) .* pr_bundl + V .* (1 .- pr_bundl);

    @assert max(abs.(V_bundl[game.S[:,1:4].==0])...) == 0
    return V_bundl, pr_bundl
end


end
