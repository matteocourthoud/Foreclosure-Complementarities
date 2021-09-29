"""Post-processing functions"""

module postprocess

include("init.jl")

using SparseArrays, LinearAlgebra, DataFrames, CSV

"""Compute transition probabilities"""
function compute_T(game)::SparseMatrixCSC{Float64,Int64}

    # Set indexes and probabilities
    idxs = (game.idx_up, game.idx_merger, game.idx_entry, game.idx_exit)
    probs = (game.Q, game.pr_merger, game.pr_entry, game.pr_exit)

    # Compute T
    T = sparse(I, game.ms, game.ms);
    for k=1:4
        K = size(idxs[k], 2);
        i = repeat(1:game.ms, outer=[K]);
        j = reshape(idxs[k][:,:,1], (game.ms*K));
        j = (j .% game.ms)
        j[j.==0] .= game.ms
        values = reshape(probs[k], (game.ms*K));
        T_temp = sparse(i, j, values, game.ms, game.ms)
        T_temp += Diagonal(vec(1 .- sum(T_temp,dims=2)))
        T = T * T_temp
    end
    # CHECK: sum(T, dims=2)
    return T
end

"""Compute and save transitions for flowchart"""
function compute_transitions(game)

    # Get transition matrix
    T = compute_T(game)

    # Fill in the summary statistics
    times = [0, 1, 3, 5, 10, 30, 100]
    states = unique(game.markets)
    block_mat = game.markets .== states'

    # Init state distribution
    distr = zeros(length(states), size(T,1))
    for (i,s) in enumerate(states)
        distr[i, s] = 1
    end

    # Init dataframe
    K = length(game.marketnames)
    colnames = ["s"; "t"; "pr_s1"; [string("q", s) for s in 1:K]]
    df = DataFrame(zeros(0, 3+K), colnames)

    # Iterate
    for i=2:length(times)
        distr *= T^(times[i] - times[i-1])
        pr_states = sum(distr, dims=1) * block_mat  ./ length(states)
        Q = distr * block_mat
        data = [states ones(length(states), 1).*times[i] pr_states' Q]
        append!(df, DataFrame(data, colnames))
    end
    # CHECK: combine(groupby(df, :t), names(df, Not(:t)) .=> sum, renamecols=false)

    # Export
    CSV.write(string("output/transitions/", game.filename, ".csv"), df, append=false, header=true)
end

"""Compute and save time series for timeline chart"""
function compute_timelines(game)

    # Get transition matrix
    T = compute_T(game)

    # Compute variables
    prices = vec(sum(game.P .* game.D, dims=2));
    profits = vec(sum(game.PI, dims=2));
    surplus = vec(game.CS);
    welfare = profits + surplus;
    V = [prices profits surplus welfare]

    # Init dataframe
    states = unique(game.markets)
    colnames = ["s", "t", "prices", "profits", "surplus", "welfare"]
    df = DataFrame(zeros(0,6), colnames)

    # Add competitive and collusive state
    smax1 = contains(game.policy, "nolearning") ? 1 : game.smax[1]
    smax2 = contains(game.policy, "nolearning") ? 1 : game.smax[2]
    omax = contains(game.policy, "nomergers") ? 0 : 1
    collusive_s = Int8.([smax1 0 smax2 0 omax])
    competitive_s = Int8.([smax1 smax1 smax2 smax2 omax])
    for s in [collusive_s, competitive_s]
        s_index = init.find_states(s, game.S)[1]
        append!(df, DataFrame([s_index 0 reshape(V[s_index,:], (1,4))], colnames))
    end

    # Init state distribution
    distr = zeros(length(states), size(T,1))
    for (i,s) in enumerate(states)
        distr[i, s] = 1
    end

    # Compute timelines
    for t=1:20
        append!(df, DataFrame([states ones(length(states)) .* t  distr * V], colnames))
        distr *= T
    end

    # Export
    CSV.write(string("output/timeseries/", game.filename, ".csv"), df, append=false, header=true)
end

"""Compute expected discounted value of variable over t periods"""
function compute_edv(v, t, s::Int64, T::SparseMatrixCSC{Float64,Int64}, beta::Float64)::Float64
    exp_values = zeros(t)
    distr_t = zeros(1, size(T,1))
    distr_t[s] = 1
    for i=1:t
        exp_values[i] = (distr_t * v)[1]
        distr_t *= T
    end
    beta = beta.^(0:t-1)
    exp_discounted_value = exp_values' * beta ./ sum(beta)
    return exp_discounted_value
end

"""Compute expected cumulative value of variable over t periods"""
function compute_ecv(v, t, s::Int64, T::SparseMatrixCSC{Float64,Int64})::Float64
    exp_values = zeros(t)
    distr_t = zeros(1, size(T,1))
    distr_t[s] = 1
    for i=1:t
        exp_values[i] = (distr_t * v)[1]
        distr_t *= T
    end
    exp_cumulative_value = sum(exp_values)
    return exp_cumulative_value
end

"""Compute expected value of variable in period t"""
function compute_ev(v, t, s::Int64, T::SparseMatrixCSC{Float64,Int64})::Float64
    distr_t = zeros(1, size(T,1))
    distr_t[s] = 1
    for i=1:t
        distr_t *= T
    end
    exp_value = (distr_t * v)[1]
    return exp_value
end

"""Get summary statistics"""
function get_sumstats(game)::DataFrame

    # Compute transitions
    T = compute_T(game);

    # Compute variables
    margin = vec(sum(game.PI[:,1:2], dims=2));
    bcost = margin .< 0;
    entry = vec(sum(game.pr_entry, dims=2));
    exit = vec(sum(game.pr_exit, dims=2));
    merger = game.S[:,5] .> 0;
    profits = vec(sum(game.PI, dims=2));
    surplus = vec(game.CS);
    welfare = profits + surplus;
    mpolyA = vec(sum(game.S[:,1:2].>0, dims=2) .== 1);
    mpolyB = vec(sum(game.S[:,3:4].>0, dims=2) .== 1);
    mpoly = vec(sum(game.S[:,1:4].>0, dims=2) .== 2);

    # Fill in the summary statistics
    states = unique(game.markets)
    beta = game.beta
    stats = zeros(length(states), 11)
    short_run = max(game.smax...)
    long_run = 1000
    for (i,s) in enumerate(states)
        stats[i,:] = [compute_edv(margin, short_run, s, T, beta),
                        compute_edv(bcost, short_run, s, T, beta),
                        compute_ecv(entry, short_run, s, T),
                        compute_ecv(exit, short_run, s, T),
                        compute_ev(merger, short_run, s, T),
                        compute_edv(profits, long_run, s, T, beta),
                        compute_edv(surplus, long_run, s, T, beta),
                        compute_edv(welfare, long_run, s, T, beta),
                        compute_ev(mpoly, long_run, s, T),
                        compute_ev(mpolyA, long_run, s, T),
                        compute_ev(mpolyB, long_run, s, T)]
    end

    # Convert to dataframe
    colnames = ["margin", "bcost", "entry", "exit", "merger", "profits", "surplus", "welfare", "mpoly", "mpolyA", "mpolyB"]
    sumstats = DataFrame(Float16.(stats), colnames)
    insertcols!(sumstats, 1, :market => game.marketnames)
    insertcols!(sumstats, 1, :sigma => game.sigma)
    insertcols!(sumstats, 1, :p0 => game.p0)
    insertcols!(sumstats, 1, :gamma => game.gamma)
    insertcols!(sumstats, 1, :c => game.c)
    insertcols!(sumstats, 1, :beta => game.beta)
    insertcols!(sumstats, 1, :alpha => game.alpha)
    return sumstats
end

"""Make all summary statistics"""
function save_to_csv(df, folder, policy)
    output_dir = string("output/", folder, "/")
    filename = string(policy, ".csv")
    file_exists = filename in readdir(output_dir)
    CSV.write(string(output_dir, filename), df, append=file_exists, header=(file_exists==false))
end

"""Generate summary statistics"""
function compute_sumstats()
    # Precompile stuff
    precompile(init.import_game, (String,))
    precompile(compute_T, (Main.init.model,))
    precompile(get_sumstats, (Main.init.model, SparseMatrixCSC{Float64,Int64}))
    precompile(save_to_csv, (Main.init.model, String, String))

    # Remove existing
    [rm("output/sumstats/$f") for f in readdir("output/sumstats/") if occursin("csv", f)]

    # Import and save
    filenames = [f for f in readdir("output/games/") if occursin("json", f)]
    for filename in filenames
        print("\n\n", filename);
        game = init.import_game(filename)
        sumstats = get_sumstats(game);
        print("\n", sumstats[:,[7:8; 12:18]], "\n\n");
        save_to_csv(sumstats, "sumstats", game.policy)
    end
end

"""Get summary statistics"""
function get_statestats(game)::DataFrame

    # Compute last state of each market
    marketstates = [sum(game.S[:,1:4].>0, dims=2) game.S[:,5]];
    states = unique([findlast(all(marketstates.==marketstates[row,:]', dims=2))[1] for row in 1:game.ms])

    # Compute variables
    margin = vec(sum(game.PI[:,1:2], dims=2));
    bcost = margin .< 0;
    entry = vec(sum(game.pr_entry, dims=2));
    exit = vec(sum(game.pr_exit, dims=2));
    merger = game.S[:,5] .> 0;
    profits = vec(sum(game.PI, dims=2));
    surplus = vec(game.CS);
    welfare = profits + surplus;
    mpoly = vec(sum(game.S[:,1:4].>0, dims=2) .== 2);

    # Compute stats
    stats = [margin bcost entry exit merger profits surplus welfare];
    stats = stats[states,:];

    # Make it a dataframe
    colnames = ["margin", "bcost", "entry", "exit", "merger", "profits", "surplus", "welfare"]
    statestats = DataFrame(Float16.(stats), colnames)
    insertcols!(statestats, 1, :market => game.marketnames)
    insertcols!(statestats, 1, :sigma => game.sigma)
    insertcols!(statestats, 1, :alpha => game.alpha)

    return statestats
end

end
