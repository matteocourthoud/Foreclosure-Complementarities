"""Functions for the learning-by-doing model"""

module model_lbd

using NLsolve, Optim

include("init.jl")
include("dynamics.jl")
include("postprocess.jl")

"""Initialize prices"""
function init_PV(game)::Tuple{Matrix{Float64}, Matrix{Float64}}
    if game.policy == "limitedmergers"
        filename = "nomergers_a$(Int(alpha*100))g$(Int(gamma*100))s$(Int(sigma*10)).json"
        temp = init.import_game(filename)
        P = temp.P
        V = game.V
        print("\nLoading alternative P and V")
    else
        P = update_P_BR(game, compute_W(game, game.V))
        V = game.V
    end
    return P, V
end

"""Correct prices for 1 firm with 1 product: equal split"""
function correct_P(game, P::Matrix{Float64})::Matrix{Float64}
    for row=1:size(game.S,1)
        s = game.S[row,:]
        mc = game.mc[row,:]
        p = P[row,:]
        # Equally split surplus in symmetric states
        if (s[6]>0)
            surplus = (p[1] + p[3]) - (mc[1] + mc[3])
            p[[1,3]] = surplus/2 .+ mc[[1,3]]
            if (s[6]>0) && (s[5]==3)
                surplus = (p[2] + p[4]) - (mc[2] + mc[4])
                p[[2,4]] = surplus/2 .+ mc[[2,4]]
            end
        end
        P[row,:] = p
    end
    return P
end

"""Precompute stuff"""
function precompute(game, row::Int, W::Array{Float64,3})

    # Check active firms
    active_n = game.active_firms[row, :]
    active_out = game.active_outcomes[row,:].>0
    out = game.active_outcomes[row, active_out.>0]'

    # Marginal cost and partner
    mc = game.mc[row, active_n[1:4]]
    partner = game.ownership[game.S[row,5]+1, active_n[1:4]]
    if (partner[1]>0) && (game.S[row, 2]==0)
        partner[1] -= 1
    end

    # Outcomes and value
    outcomes = game.outcomes[active_out, active_n]'
    w = W[row, active_out, active_n[1:4]]'
    sign = (outcomes[1:end-1,:] .> 0) .- (outcomes[1:end-1,:] .== 0)
    w_signed = copy(w)
    w_signed += W[row, active_out, game.partner[active_n[1:4]]]' .* (partner.>0)
    w_signed .*= sign
    joint_p = Int.([1,2,1,2][active_n[1:4]])

    return active_n, out, mc, outcomes, partner, w, w_signed, joint_p
end

"""Compute demand"""
function demand(p::Vector{Float64}, sigma::Float64, p0::Float64, outcomes, out)::Tuple{Matrix{Float64},Vector{Float64}}
    u = - [p; p0] .* sigma          # Utility
    u = u .- max(u...)              # Normalize
    e = exp.(u' * outcomes) .* out  # Exponential of utility * outcomes
    q = e ./ sum(e)                 # Demand of each system (should sum to 1)
    d = outcomes * q'               # Demand of each product (should NOT sum to 1)
    return q, d[1:end-1]
end

"""Compute profits for BR iteration"""
function BR(x::Vector{Float64}, n::Int, p::Vector{Float64}, sigma::Float64, p0::Float64, mc::Vector{Float64}, outcomes, partner::Vector{Int}, w, out)::Float64
    p[n] = x[1]                                         # Insert price of firm n
    q, d = demand(p, sigma, p0, outcomes, out)# Compute demand
    V = (p .- mc) .* d + sum(w .* q, dims=2)            # Compute value
    V[n] += (partner[n] > 0) ? V[partner[n]] : 0        # Correct for ownership
    return V[n]
end

"""Compute price by best reply itearation"""
function update_p_BR(game, p::Vector{Float64}, row::Int, W::Array{Float64,3})::Vector{Float64}

    # Precompute things
    active_n, out, mc, outcomes, partner, w, w_signed, joint_p = precompute(game, row, W)

    # Iterate best reply
    dist = 1
    iter = 0
    br = copy(p[active_n[1:4]])
    while (dist > game.accuracy) && (iter<100)
        for n=1:length(br)
            br[n] = optimize((x -> -BR(x, n, br, game.sigma, game.p0, mc, outcomes, partner, w, out)), [br[n]], LBFGS()).minimizer[1]
        end
        dist = max(abs.(br .- p[active_n[1:4]])...)
        p[active_n[1:4]] = copy(br);
        iter += 1
    end
    return p
end

"""Update prices in all states"""
function update_P_BR(game, W::Array{Float64,3})::Matrix{Float64}
    P = copy(game.P)
    for row=1:size(game.S,1)
        p = copy(game.P[row, :])
        P[row, :] = update_p_BR(game, p, row, W)
    end
    P = correct_P(game, P)
    return P
end

"""First order condition"""
function FOC(p::Vector{Float64}, sigma::Float64, p0::Float64, mc::Vector{Float64}, outcomes, partner::Vector{Int}, w_signed::Matrix{Float64}, joint_p::Vector{Int}, out)::Vector{Float64}

    # Compute demand
    q, d = demand(p, sigma, p0, outcomes, out)

    # Compute extra: future value
    dd = ((1 .- d) .* (outcomes[1:end-1,:] .> 0) .+ d .* (outcomes[1:end-1,:] .== 0) )
    EV = vec(sum(dd .* q .* w_signed, dims=2))

    # Compute extra: ownership
    EO = zeros(length(p))
    for (n, ptn) in enumerate(partner)
        if ptn>0
            EO[n] = (q[joint_p[n]] - d[n]*d[ptn]) * (p[ptn] .- mc[ptn])
        end
    end

    # Compute zero
    z = d .* (1 .- d) .* (p .- mc) .+ EV .+ EO .- d ./ sigma
    return z
end

"""Update p by numerically solving first order condition"""
function update_p_FOC(game, W::Array{Float64,3}, row::Int)::Vector{Float64}

    # Precompute things
    active_n, out, mc, outcomes, partner, w, w_signed, joint_p = precompute(game, row, W)

    # Init
    p = game.P[row, :]

    # Solve
    solution = nlsolve((x -> FOC(x, game.sigma, game.p0, mc, outcomes, partner, w_signed, joint_p, out)), p[active_n[1:4]])
    p[active_n[1:4]] = solution.zero

    # Compute demand
    q, d = demand(p[active_n[1:4]], game.sigma, game.p0, outcomes, out)

    # Consider model not solved if not zero, there are nan prices or zero demand
    threshold = max(1-game.sigma/10, 0.3)
    not_solved = (!solution.f_converged) || (q[end]>threshold) #|| (max(p...)>2*game.p0)
    if not_solved
        p = update_p_BR(game, 2 .* game.mc[row,:], row, W)
    end

    return p
end

"""Update prices in all states"""
function update_P(game, W::Array{Float64,3})::Matrix{Float64}
    P = copy(game.P)
    for row=1:size(game.S,1)
        P[row,:] = update_p_FOC(game, W, row)
    end
    P = correct_P(game, P)
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
    U = -[P ones(game.ms,1).*game.p0] .* game.sigma     # Consumer utility
    E = exp.(U * game.outcomes') .* game.active_outcomes# Expontential utility
    Q = E ./ (sum(E, dims=2))                           # Product demand
    D = Q * game.outcomes[:,1:4]                        # Firm demand
    PI = D .* (P .- game.mc)                            # Profits
    CS = 2*game.p0 .+ log.(sum(E, dims=2)) ./ game.sigma
    return Q, D, PI, CS
end

"""Update value function"""
function update_V_pricing(game, V::Matrix{Float64}, args...)
    W = (length(args)==1) ? compute_W(game, args[1]) : compute_W(game, V)
    P = (length(args)==2) ? args[2] : update_P(game, W)
    # Update values
    Q, D, PI, CS = compute_PI(game, P)
    V1 = Float64.(zeros(size(V)))
    for n=1:4
        V1[:,n] = PI[:,n] + sum(Q .* W[:,:,n], dims=2);
    end
    return V1, P, Q, D, PI, CS
end

end
