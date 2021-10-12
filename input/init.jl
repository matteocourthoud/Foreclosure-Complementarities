"""Model of algorithms and competition"""

module init

using Statistics, StructTypes, Optim, JSON3

"""Model properties"""
Base.@kwdef mutable struct model

    """Default Properties"""
    policy::String = "baseline"                 # Name of policy
    smax::Vector{Int8} = Int8[5,1]              # Number of states per side of the market
    alpha::Float64 = 0.8                        # Learning parameter
    beta::Float64 = 0.95                        # Discount factor
    c::Float64 = 1.0                            # Marginal cost
    gamma::Float64 = 0.0                        # Complementarity parameter
    sigma::Float64 = 5.0                        # Competition parameter
    p0::Float64 = 1.5                           # Value of data
    accuracy::Float64 = 1e-8                    # Approximation accuracy
    cost_entry::Vector{Float64} = [0;10];       # Entry cost
    value_exit::Vector{Float64} = [0;0.5];      # Exit scrap values
    cost_merger::Vector{Float64} = [0;10];      # Merger cost
    entry::Bool = true                          # Entry
    exit::Bool = true                           # Exit
    mergers::Bool = (policy != "nomergers")     # Mergers
    filename::String = string(policy, "_a", Int64(floor(alpha*100)), "g", Int64(floor(gamma*100)), "s", Int64(floor(sigma*10)))
    verbose::Bool = true

    """Precomputed stuff to speed up computation"""
    firms::Vector{Int8} = Int8[1, 2, 3, 4]      # Firms, sigma=5.0
    rival::Vector{Int8} = Int8[2, 1, 4, 3]      # Rivals
    partner::Vector{Int8} = Int8[3, 4, 1, 2]    # Partners
    merger_pairs::Matrix{Int8} = Int8[1 3; 2 4; 1 4; 2 3]; # Merger pairs
    outcomes::Matrix{Int8} = Int8[1 0 1 0 0
                                  0 1 0 1 0
                                  1 0 0 1 0
                                  0 1 1 0 0
                                  1 0 0 0 1
                                  0 1 0 0 1
                                  0 0 1 0 1
                                  0 0 0 1 1
                                  0 0 0 0 2]    # Sale outcomes

    ownership::Matrix{Int8} = Int8[0 0 0 0
                                   3 0 1 0
                                   0 4 0 2
                                   3 4 1 2]     # Ownership matrix

    """Derived Properties"""
    S::Matrix{Int8} = get_state_space(smax, policy) # State space
    ms::Int64 = size(S,1)                       # Dimension of the state space
    mc::Matrix{Float64} = compute_mc(S, alpha, c, smax, policy) # Marginal cost
    P::Matrix{Float64} = 2 .* mc                # Prices
    V::Matrix{Float64} = zeros(ms,4)            # Value function
    D::Matrix{Float64} = zeros(ms,4)            # Demand (per firm)
    Q::Matrix{Float64} = zeros(ms,size(outcomes,1))  # Demand (per system)
    PI::Matrix{Float64} = zeros(ms,4)           # Profits
    CS::Matrix{Float64} = zeros(ms,1)           # Consumer surplus
    pr_entry::Matrix{Float64} = zeros(ms,4)     # Probability of entry
    pr_exit::Matrix{Float64} = zeros(ms,4)      # Probability of exit
    pr_merger::Matrix{Float64} = zeros(ms,4)    # Probability of merger
    active_firms::BitMatrix = Bool.([S[:,1:4] .> 0 ones(ms, 1)])
    active_outcomes::Matrix{Float64} = compute_active_outcomes(S, ms, outcomes, gamma, policy)
    markets::Vector{Int64} = compute_markets(S, ms)
    marketnames::Vector{String} = compute_marketnames(S)
    idx_up::Array{Int64,3} = compute_idx_up(S, ms, smax, outcomes, active_outcomes, policy)
    idx_entry::Array{Int64,3} = compute_idx_entry(S, ms, entry)
    idx_exit::Array{Int64,3} = compute_idx_exit(S, ms, exit, rival, ownership)
    idx_merger::Array{Int64,3} = compute_idx_merger(S, ms, mergers, rival, merger_pairs, policy)

end

# Save struct type
StructTypes.StructType(::Type{model}) = StructTypes.Mutable()

"""Get actions of the platform: order of firms"""
function get_state_space(smax::Vector{Int8}, policy::String)::Matrix{Int8}
    S = Int8.(zeros(0,5))
    i1max = (policy == "nolearning") ? 1 : smax[1]
    i3max = (policy == "nolearning") ? 1 : smax[2]
    O = (policy == "nomergers") ? [0] : [0,1,3]
    for o=O
        for i1=1:i1max
            for i3=1:i3max
                i2max = i1 * (o!=1) + i1max * (o==1)
                for i2=(o==3):i2max
                    i4max = i3 * (o==0) + i3max * (o!=0)
                    for i4=(o==3):i4max
                        s = reshape([i1, i2, i3, i4, o], (1,5))
                        S = [S; s]
                    end
                end
            end
        end
    end
    # Sort columns
    S = [S reshape(sum(S[:,1:4].>0, dims=2), (size(S, 1), 1))]
    S = sortslices(S, dims=1, by=x->x[3], rev=false)
    S = sortslices(S, dims=1, by=x->x[1], rev=false)
    S = sortslices(S, dims=1, by=x->x[2], rev=false)
    S = sortslices(S, dims=1, by=x->x[4], rev=false)
    S = sortslices(S, dims=1, by=x->x[5], rev=false)
    S = sortslices(S, dims=1, by=x->x[6], rev=false)
    S = S[:, 1:5]
    return S
end

"""Find index of a single state"""
function find_state(s::Vector{Int8}, S::Matrix{Int8})::Vector{Int64}
    row_index = findall(all(reshape(s, (1,5)) .== S, dims=2))[1][1]
    row_indexes = LinearIndices(S)[row_index, 1:4]
    return row_indexes
end

"""Find index of a set of states"""
function find_states(states::Matrix{Int8}, S::Matrix{Int8})::Matrix{Int64}

    # Init
    indexes = Int64.(zeros(size(states,1),4));

    # Loop over states
    for i=1:size(states,1)

        s = states[i,:];

        # Try to find the state
        order = [1,2,3,4];
        j = findfirst(all(reshape(s, (1,5)) .== S, dims=2));

        # If not found, change order
        if isnothing(j)
            # If 2/4 ownership, swap
            if s[5]==2
                order = [2,1,4,3];
                s[5] = 1;
            # If lower states and no ownership, double swap
            elseif (s[1]<s[2]) && (s[3]<s[4]) && (s[5]==0)
                order = [2,1,4,3];
            # If lower states and doble ownership, double swap
            elseif ((s[1]<s[2]) || (s[3]<s[4])) && (s[5]==3)
                order = [2,1,4,3];
            # If state 3 lower and no ownership, swap
            elseif (s[3]<s[4]) && (s[5]==0)
                order = [1,2,4,3];
            # If state 1 lower and no ownership, swap
            elseif (s[1]<s[2]) && (s[5]==0)
                order = [2,1,3,4];
            end

            # Change order of state
            s = [s[order]; s[5]]
        end

        # Find index
        idx = find_state(s, S)

        # Switch back
        indexes[i,:] = idx[order]
    end
    return indexes
end

"""Compute marginal cost"""
function compute_mc(S::Matrix{Int8}, alpha::Float64, c::Float64, smax::Vector{Int8}, policy::String)::Matrix{Float64}
    mc = c .* S[:,1:4].^log2(alpha);
    if (policy == "nolearning")
        mc[:,1:2] =  c .* (S[:,1:2] .* smax[1]).^log2(alpha);
        mc[:,3:4] =  c .* (S[:,3:4] .* smax[2]).^log2(alpha);
    end
    mc[S[:,1:4].==0] .= 0
    return mc
end

"""Compute index when making a sale"""
function compute_idx_up(S::Matrix{Int8}, ms::Int64, smax::Vector{Int8}, outcomes::Matrix{Int8}, active_outcomes::Matrix{Float64}, policy::String)::Array{Int64, 3}

    # Init I_up: state x product x firm
    K = size(outcomes,1)
    idx_up = Int64.(ones(ms, K, 4))

    # Init maximum achievable state
    Sfirms = S[:,1:4];
    Smax = (policy == "nolearning") ? [1 1 1 1] : [smax[1] smax[1] smax[2] smax[2]]

    # Loop over outcomes
    for k=1:K
        # Generate next state
        out = outcomes[k,1:4];
        active_out = active_outcomes[:,k] .> 0;
        up = (Sfirms.<Smax) .* reshape(out, (1,4)) .* active_out;
        S2 = Int8.([Sfirms + up S[:, 5]])

        # Data sharing: follower is never more than one step behind
        if (policy == "datasharing")
            lag = max.(S2[:,[2,1,4,3]] - S2[:,1:4] .- 1, 0) .* (Sfirms.>0);
            S2[:,1:4] = S2[:,1:4] + lag;
        end

        # Find next state
        idx = find_states(S2, S);

        # Add state index
        for n=1:4
            idx_up[:,k,n] = idx[:,n];
        end
    end
    return idx_up
end

"""Compute entry index: in each state, to which state would each
   firm move to, for each possible firm entry?"""
function compute_idx_entry(S::Matrix{Int8}, ms::Int64, entry::Bool)::Array{Int64,3}

    # Init map. Dimensions: states x entrant x firms
    idx_entry = Int64.(zeros(ms, 4, 4));

    # Loop over entrants
    for e=1:4
        # Loop over rows
        for row=1:ms

            # Init
            s = reshape(S[row,:], (1,5));

            # Check if entry is possible
            entry_possible = (s[e]==0) && entry;

            # Update state if entry is possible
            if entry_possible
                s[e] = 1;
            end

            # Find state
            idx_entry[row,e,:] = init.find_states(s, S);
        end
    end
    # CHECK
    @assert min(idx_entry...) > 0
    return idx_entry
end

"""Compute exit index: in each state, to which state would each
   firm move to, for each possible firm entry?"""
function compute_idx_exit(S::Matrix{Int8}, ms::Int64, exit::Bool, rival::Vector{Int8}, ownership::Matrix{Int8})::Array{Int64,3}

    # Init map. Dimensions: states x exiters x firms
    idx_exit = Int64.(zeros(ms, 4, 4));

    # Loop over exiters
    for e=1:4
        # Loop over states
        for row=1:ms

            # Init
            s = reshape(S[row,:], (1,5));
            o = s[5];
            partner = ownership[o+1, e]

            # Check if exit is possible
            # - both firm and its rival are active
            # - if firm has partner, there must be 4 firms
            exit_possible = (min([s[e], s[rival[e]]]...)>0) && exit;
            if partner>0
                exit_possible = exit_possible && (sum(s[1:4].>0)==4)
            end

            # If exit possible
            if exit_possible
                # Set exiter state equal to zero
                s[e] = 0;
                # If firm has partner, they both exit
                if partner>0
                    s[partner] = 0;
                end

                # Decrease ownership if nonzero
                s[5] = s[5] - (e in [1,3]) * (o in [1,3]);
                s[5] = s[5] - 2 * (e in [2,4]) * (o == 3);

            end

            # Find state
            idx_exit[row,e,:] = init.find_states(s, S);
        end
    end
    # CHECK
    @assert min(idx_exit...) > 0
    return idx_exit
end

"""Compute merger index: in each state, to which state would each
   firm move to, for each possible merger?"""
function compute_idx_merger(S::Matrix{Int8}, ms::Int64, mergers::Bool, rival::Vector{Int8}, merger_pairs::Matrix{Int8}, policy::String)::Array{Int64,3}

    # Init map. Dimensions: states x merger pairs x firms
    idx_merger = Int64.(zeros(ms, 4, 4));

    # Loop over merger pairs
    for p=1:4
        # Loop over states
        for row=1:ms

            # Init
            s = reshape(S[row,:], (1,5));
            o = s[5];

            # Check if merger possible
            # - all merging parties active
            # - either no ownership or ownership and 2nd pair (2-4)
            merging_firms = merger_pairs[p,:]
            merger_possible = min(s[merging_firms]...)>0;
            merger_possible = merger_possible && ((o==0) || ((o==1) && (p==2)));
            merger_possible = merger_possible && mergers;

            # Merger policy: no merger between market leaders
            if (policy == "limitedmergers")
                rival_firms = rival[merger_pairs[p,:]]
                merger_possible = merger_possible && all(s[merging_firms] .<= s[rival_firms]);
            end

            # If merger is possible
            if merger_possible

                # If merger pair is 3/4 and no own, swap 3-4
                if (p in [3,4]) && (o==0) && merger_possible
                    s[1:4] = s[[1,2,4,3]];
                end

                # Change ownership
                s[5] = s[5] + (p in [1,3]) * (o in [0,2]);
                s[5] = s[5] + 2 * (p in [2,4]) * (o in [0,1]);
            end

            # Find state
            idx_merger[row,p,:] = init.find_states(s, S);

            # If merger pair is 3/4 and no own, swap back idx 3-4
            if (p in [3,4]) && (o==0) && merger_possible
                idx_merger[row,p,1:4] = idx_merger[row,p,[1,2,4,3]];
            end
        end
    end
    # CHECK
    @assert min(idx_merger...) > 0
    return idx_merger
end

"""Compute which outcomes are active in each state"""
function compute_active_outcomes(S::Matrix{Int8}, ms::Int64, outcomes::Matrix{Int8}, gamma::Float64, policy::String)::Matrix{Float64}
    active_outcomes = zeros(ms, size(outcomes,1))
    # Both firms must be active
    for i=1:ms
        active_outcomes[i,:] = (sum([S[i, 1:4].>0; 1] .* outcomes', dims=1).==2)
    end
    # No bundling: no limits to products combinations
    if policy != "nobundling"
        # With double ownership, only 1,2,9 are active
        active_outcomes[S[:,5] .== 3, 3:end-1] .= 0
        # With single ownership only 1,2,6,8,9 are active
        if policy != "limitedbundling"
            active_outcomes[S[:,5] .== 1, [3,4,5,7]] .= 0
        end
    end
    # Partial complementarity
    active_outcomes[:, 5:end-1] = active_outcomes[:, 5:end-1] .* gamma
    return active_outcomes
end

"""Compute markets"""
function compute_markets(S::Matrix{Int8}, ms::Int64)::Vector{Int64}
    marketstates = [sum(S[:,1:4].>0, dims=2) S[:,5]];
    markets = [findfirst(all(marketstates.==marketstates[row,:]', dims=2))[1] for row in 1:ms]
    return markets
end

"""Compute market names"""
function compute_marketnames(S::Matrix{Int8})::Vector{String}
    marketstates = [sum(S[:,1:4].>0, dims=2) S[:,5]];
    uniquemarkets = unique(marketstates, dims=1)
    marketnames = [join(string.(uniquemarkets[row,:])) for row in 1:size(uniquemarkets,1)]
    return marketnames
end

"""Import game from json file"""
function import_game(filename::String)
    # Read the json
    path = string("output/games/", filename)
    json_obj = JSON3.read(read(path, String))
    # Initialize game
    policy = json_obj.policy
    game = model(policy=policy)
    # Reshape attributes
    for field in fieldnames(typeof(game))
        value = json_obj[field]
        if typeof(value) != String
            fieldsize = size(getfield(game, field))
            if length(fieldsize)>0
                value = reshape(value, fieldsize)
            end
        end
        setproperty!(game, field, value)
    end
    return game
end

end
