"""Network Effects Functions"""

module lbd

using NLsolve, Optim, JSON3

include("init.jl")
include("dynamics.jl")
include("postprocess.jl")

"""Export game"""
function export_game(game)
    filename = game.filename
    open("output/games/$filename.json", "w") do io
        JSON3.write(io, game)
    end
end

"""Initialize price"""
function init_P(game)::Tuple{Matrix{Float64},Matrix{Float64}}
    game.P = - 0.5 * game.scale
    V = game.scale
    W = compute_W(game, V)
    P = update_P_BR(game, W)
    return P, V

end

"""Correct prices for 1 firm with 1 product"""
function correct_P(game, P::Matrix{Float64})::Matrix{Float64}
    for row=1:size(game.S,1)
        s = game.S[row,:]
        scale = game.scale[row,:]
        p = P[row,:]
        sym_state = ((sum(s[1:4] .> 0).==2) && (s[5]==1)) || (s[5]==3);
        if sym_state
            surplus = (scale[1] + scale[3]) + (p[1] + p[3])
            p[[1,3]] = surplus/2 .- scale[[1,3]]
            if (s[5]==3)
                surplus = (scale[2] + scale[4]) + (p[2] + p[4])
                p[[2,4]] = surplus/2 .- scale[[2,4]]
            end
        end
        P[row,:] = p
    end
    return P
end

"""Compute demand"""
function demand(p::Vector{Float64}, sigma::Float64, outcomes, out)::Tuple{Matrix{Float64},Vector{Float64}}
    u = - [p; 0] .* sigma       # Utility
    u = u .- max(u...)          # Normalize
    e = exp.(u' * outcomes)     # Exponential of utility
    q = e ./ sum(e)             # Demand of each system (should sum to 1)
    d = outcomes * q'           # Demand of each product (should NOT sum to 1)
    return q, d[1:end-1]
end

"""Compute profits for BR iteration"""
function BR(x::Vector{Float64}, n::Int64, p::Vector{Float64}, game, scale::Vector{Float64}, outcomes, partner::Vector{Int8}, w, out)::Float64

    # Insert price of firm n
    p[n] = x[1]

    # Compute demand
    q, d = demand(p, game.sigma, outcomes, out)

    # Compute value
    V = (p .+ scale) .* d + sum(w .* q, dims=2)

    # Compute value of firm n
    Vn = V[n]

    # Correct for ownership
    if partner[n] > 0
        Vn += V[partner[n]]
    end
    return Vn
end

"""Compute price by best reply itearation"""
function update_p_BR(game, p, row, W)::Vector{Float64}

    # Check active firms
    active_n = game.active_firms[row,:]
    active_out = game.active_outcomes[row,:].>0
    out = game.active_outcomes[row, active_out.>0]'

    # Init
    scale = game.scale[row, active_n[1:4]]
    partner = game.ownership[game.S[row,5]+1, active_n[1:4]]
    if (partner[1]>0) & (game.S[row, 2]==0)
        partner[1] -= 1
    end
    outcomes = game.outcomes[active_out, active_n]'
    w = W[row, active_out, active_n[1:4]]'

    # Iterate best reply
    dist = 1
    iter = 0
    br = copy(p[active_n[1:4]])
    while (dist > game.accuracy) && (iter<100)
        for n=1:length(br)
            br[n] = optimize((x -> -BR(x, n, br, game, scale, outcomes, partner, w, out)), [br[n]], LBFGS()).minimizer[1]
        end
        dist = max(abs.(br .- p[active_n[1:4]])...)
        p[active_n[1:4]] = copy(br);
        iter += 1
    end
    return p
end

"""Update prices in all states"""
function update_P_BR(game, W)::Matrix{Float64}
    P = zeros(size(game.P))
    for row=1:size(game.S,1)
        p = copy(game.P[row, :])
        P[row, :] = update_p_BR(game, p, row, W)
    end
    return P
end

"""First order condition"""
function FOC(p::Vector{Float64}, game, scale::Vector{Float64}, outcomes, partner::Vector{Int8}, w_signed, joint_p::Vector{Int8}, out)::Vector{Float64}

    # Compute demand
    q, d = demand(p, game.sigma, outcomes, out)

    # Compute extra: future value
    dd = ((1 .- d) .* (outcomes[1:end-1,:] .> 0) .+ d .* (outcomes[1:end-1,:] .== 0) )
    EV = vec(sum(dd .* q .* w_signed, dims=2))

    # Compute extra: ownership
    EO = zeros(length(p))
    for (n, n_) in enumerate(partner)
        if n_>0
            EO[n] = (q[joint_p[n]] - d[n]*d[n_]) * (p[n_] .+ scale[n_])
        end
    end

    # Compute zero
    z = d .* (1 .- d) .* (p .+ scale) .+ EV .+ EO .- d ./ game.sigma
    return z
end

"""Update p by numerically solving first order condition"""
function update_p_FOC(game, W, row)::Vector{Float64}

    # Check active firms
    active_n = game.active_firms[row, :]
    active_out = game.active_outcomes[row,:].>0
    out = game.active_outcomes[row, active_out.>0]'

    # Init
    scale = game.scale[row, active_n[1:4]]
    partner = game.ownership[game.S[row,5]+1, active_n[1:4]]
    if (partner[1]>0) & (game.S[row, 2]==0)
        partner[1] -= 1
    end
    outcomes = game.outcomes[active_out, active_n]'
    w = W[row, active_out, active_n[1:4]]'
    sign = (outcomes[1:end-1,:] .> 0) .- (outcomes[1:end-1,:] .== 0)
    w += (W[row, active_out, game.partner[active_n[1:4]]]' .* (partner.>0))
    w_signed = Float64.(sign .* w)
    joint_p = Int8.([1,2,1,2][active_n[1:4]])

    # Init
    p = game.P[row, :]

    # Solve
    solution = nlsolve((x -> FOC(x, game, scale, outcomes, partner, w_signed, joint_p, out)), p[active_n[1:4]])
    p[active_n[1:4]] = solution.zero

    # Compute demand
    q, d = demand(p[active_n[1:4]], game.sigma, outcomes, out)

    # Consider model not solved if not zero, there are nan prices or zero demand
    not_solved = (solution.f_converged==false) || (sum(d)<0.5) || (max(p...)>10)
    if not_solved
        p = 0 * game.scale[row,:]
        p = update_p_BR(game, p, row, W)
    end

    return p
end

"""Update prices in all states"""
function update_P(game, W::Array{Float64,3})::Matrix{Float64}
    P = copy(game.P)
    for row=1:size(game.S,1)
        P[row,:] = update_p_FOC(game, W, row)
    end
    return P
end

"""Compute outcome-specific future values"""
function compute_W(game, V::Matrix{Float64})::Array{Float64,3}
    W = zeros(game.ms, size(game.outcomes,1), 4);
    EV = game.beta * V;
    for n=1:4
        W[:,:,n] = EV[game.idx_up[:,:,n]];
    end
    return W
end

"""Compute demand and profits"""
function compute_PI(game, P::Matrix{Float64})::Tuple{Matrix{Float64},Matrix{Float64},Matrix{Float64},Array{Float64}}
    U = -[P zeros(game.ms,1)] .* game.sigma             # Consumer utility
    E = exp.(U * game.outcomes') .* game.active_outcomes# Expontential utility
    Q = E ./ (sum(E, dims=2))                           # Product demand
    D = Q * game.outcomes[:,1:4]                        # Firm demand
    PI = D .* (P .+ game.scale)                         # Profits
    bonus = 1 .+ (game.S[:,5]./3)
    CS = log.(sum(E, dims=2) .* bonus) ./ game.sigma
    return Q, D, PI, CS
end

"""Update value function"""
function update_V(game, V::Matrix{Float64}, args...)
    W = (length(args)==1) ? compute_W(game, args[1]) : compute_W(game, V)
    P = (length(args)==2) ? args[2] : update_P(game, W)
    P = correct_P(game, P) # TODO: check correction
    # Update values
    Q, D, PI, CS = compute_PI(game, P)
    V1 = Float64.(zeros(size(V)))
    for n=1:4
        V1[:,n] = PI[:,n] + sum(Q .* W[:,:,n], dims=2);
    end
    return V1, P, Q, D, PI, CS
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

    # Init values
    game.P, game.V = init_P(game)

    # Solve game
    print("\n\nSolving ", game.filename, "\n----------------------\n")
    rate = 1
    dist = 1
    iter = 0
    while (dist>100*game.accuracy) && (iter<1000)
        V4, game.pr_exit = dynamics.update_V_exit(game, game.V);
        V3, game.pr_entry = dynamics.update_V_entry(game, V4);
        V2, game.pr_merger = dynamics.update_V_merger(game, V3);
        V1, game.P, game.Q, game.D, game.PI, game.CS = update_V(game, V2);
        # Compute distance
        rate = 0.9*rate + 0.1*max(abs.(game.V - V1)...)/dist
        dist = max(abs.(game.V - V1)...);
        # Update value function
        r = rand() * (rate > 1) * (iter > 300);
        game.V = game.V .* r + V1 .* (1-r)
        # Print update
        if iter % 10 == 0
            print("\rIter ", iter, ": ", dist, " (", rate, ")")
        end
        iter += 1
    end

    # Print game
    if game.verbose
        sumstats = postprocess.get_sumstats(game);
        print("\n\n", sumstats[:,[3,4,6,7,9,10,11,12]], "\n\n");
        #statestats = postprocess.get_statestats(game);
        #print("\n\n", statestats[:,[3,4,6,7,9,10,11]], "\n\n");
    end

    # Export game
    export_game(game)

    # If there were issues, save it
    if (dist>100*game.accuracy) || (min(sum(game.D, dims=2)...)<0.01)
        open("issues.txt","a") do io
           println(io, string(game.filename, " : ", dist))
        end
    end
    return game
end

end
