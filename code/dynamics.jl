"""Entry, exit, mergers and bundling functions"""

module dynamics

"""Compute probability"""
function compute_p(Delta::Vector{Float64}, c::Vector{Float64})::Tuple{Vector{Float64},Vector{Float64}}
    # Uniform
    p = min.(max.( (Delta .- c[1]) ./ (c[2]-c[1]), 0), 1) ./ 4;
    ev = max.(c[1] .+ (min.(Delta, c[2]) .- c[1]) ./ 2, 0);
    # Exponential
    #Delta = max.(Delta, 0)
    #p = (1 .- exp.( - Delta ./ c[2])) ./ 4;
    #ev = c[2] .- Delta .* exp.( - Delta ./ c[2]) ./ (1 .- exp.( - Delta ./ c[2]))
    #ev[Delta.==0] .= 0
    return p, ev
end

"""Update value for exit"""
function update_V_exit(game, V::Matrix{Float64}, args...)::Tuple{Matrix{Float64},Matrix{Float64}}

    # Init
    V_exit = zeros(size(V));
    pr_exit = zeros(size(game.pr_exit));

    # Return if no exit
    if !game.exit
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
        Delta[rows] = (Delta[rows] + Delta_partner)/2; # TODO: double exit = double cost

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

        # Update value of firms that don't exit
        V1[:,e] = V[:,e];
        V1[rows,game.partner[e]] = V[rows,game.partner[e]];
        V_exit += V1 .* pr_exit[:,e];
    end

    # Update value with residual probability
    V_exit += V .* (1 .- sum(pr_exit, dims=2));

    # Checks
    @assert min(pr_exit...)>=-game.accuracy
    @assert max(abs.(V_exit[game.S[:,1:4].==0])...) == min(abs.(V_exit[game.S[:,1:4].==0])...) == 0
    return V_exit, pr_exit
end

"""Update value for entry"""
function update_V_entry(game, V::Matrix{Float64}, args...)::Tuple{Matrix{Float64},Matrix{Float64}}

    # Init value
    V_entry = zeros(size(V));
    pr_entry = zeros(size(game.pr_entry));

    # Return if no entry
    if !game.entry
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

    # Checks
    @assert min(pr_entry...)>=-game.accuracy
    @assert max(abs.(V_entry[game.S[:,1:4].==0])...) == min(abs.(V_entry[game.S[:,1:4].==0])...) == 0
    return V_entry, pr_entry
end

"""Update value for mergers"""
function update_V_merger(game, V::Matrix{Float64}, args...)::Tuple{Matrix{Float64}, Matrix{Float64}}

    # Init
    V_merger = zeros(size(V));
    pr_merger = zeros(size(game.pr_merger));

    # Return if no mergers
    if game.mergers == false
        return V, pr_merger
    end

    # Loop over all possible merger pairs
    for p=1:size(game.merger_pairs,1)

        # Compute merger value
        V1 = V[game.idx_merger[:,p,:]];

        # Compute policy if not provided
        c = game.cost_merger;
        Delta = reshape(sum(V1[:,game.merger_pairs[p,:]] - V[:,game.merger_pairs[p,:]], dims=2), game.ms);
        pr_merger[:,p], exp_cost = compute_p(Delta, c);

        # Update value with outcome of Nash bargaining: CHECK
        V1[:,game.merger_pairs[p,:]] = V[:,game.merger_pairs[p,:]] .+ (Delta-exp_cost) ./ 2
        V_merger += V1 .* pr_merger[:,p];
    end

    # Add residual value
    V_merger += V .* (1 .- sum(pr_merger, dims=2));

    # Checks
    @assert min(pr_merger...)>=-game.accuracy
    @assert max(abs.(V_merger[game.S[:,1:4].==0])...) == min(abs.(V_merger[game.S[:,1:4].==0])...) == 0
    return V_merger, pr_merger
end

"""Update value for bundling"""
function update_V_bundling(game, V::Matrix{Float64}, args...)::Tuple{Matrix{Float64},Matrix{Float64}}

    # Init
    V_bundling = zeros(size(V));
    pr_bundling= zeros(size(game.pr_bundling));

    # Return if no bundling
    if game.bundling == false
        return V, pr_bundling
    end

    # Update with incentives
    V_incentives = (length(args)==1) ? args[1] : V

    # Loop over all possible merger pairs
    for b=1:4

        # Compute bundling value
        V1 = V[game.idx_bundling[:,b,:]];
        V1_incentives = V_incentives[game.idx_bundling[:,b,:]];

        # Compute policy if not provided
        c = game.cost_bundling;
        pair = [b, game.partner[b]]
        Delta = reshape(sum(V1[:,pair] - V[:,pair], dims=2), game.ms);
        Delta_incentives = reshape(sum(V1_incentives[:,pair] - V_incentives[:,pair], dims=2), game.ms);
        if length(args) < 2
            pr_bundling[:,b], exp_cost = compute_p(Delta_incentives, c);
        else # CHECK but should be correct
            pr_bundling[:,b] = args[2][:,b]
            exp_cost = (c[1] .+ (c[2] .* pr_bundling[:,b] + c[1] .* (1 .- pr_bundling[:,b])) ./2 ) .* pr_bundling[:,b]
        end

        # Update value with outcome of Nash bargaining: CHECK
        V1[:,pair] = V[:,pair] .+ (Delta-exp_cost) ./ 2
        V_bundling += V1 .* pr_bundling[:,b];
    end

    # Add residual value
    V_bundling += V .* (1 .- sum(pr_bundling, dims=2));

    # Checks
    @assert min(pr_bundling...)>=-game.accuracy
    @assert max(abs.(V_bundling[game.S[:,1:4].==0])...) == min(abs.(V_bundling[game.S[:,1:4].==0])...) == 0
    return V_bundling, pr_bundling
end

end
