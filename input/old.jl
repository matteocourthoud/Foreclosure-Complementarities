module old


function BR_simple(game, A::Matrix{Float32}, B::Matrix{Float32}, n::Int8)::Float32
    """
    Compute best reply of player n
    -------------
        # DEMAND    d := A * exp(-s*x) / ( A * exp(-s*x) + B )

        # PROFITS   p := d * (x - c)

        # FOC       0 == s*B*(x - c) - (e + B)
            where   e := A * exp(-s*x)

        # BR        br = ( # + lambertw( exp(-#) * A / B) ) / s
            where   # := 1 + s*c
    """

    # Compute best reply
    arg = 1 + game.sigma*game.c
    br = ( arg + lambertw( exp(-arg) * A[n] / B[n] )) / game.sigma
    return br
end

function update_p_BR_simple(game, s::Vector{Int8}, w, p0::Vector{Float32})::Vector{Float32}
    """Update price in one state by BR iteration"""

    # Check active firms
    active = Int8.(s[1:4] .> 0)
    active_n = Int8.([n for n in 1:4 if active[n]>0])

    # Init x
    x = p0 .*  active[1:4]
    br = Float32.(zeros(4))

    # Compute exp(state) of each product
    Es = exp.(game.alpha .* s[1:4]) .* active

    # Best reply once
    for n in active_n
        Ep = exp.(-game.sigma .* x) .* active
        A = Es .* sum(((Es .* Ep)[game.BR_index[:,1:2]]), dims=2)
        p = prod([Es .* Ep; exp(2*game.sigma*game.v0)][game.outcomes_index], dims=1)
        B = sum(p[game.nonoutcomes_index], dims=2) .* active
        x[n] = BR_simple(game, A, B, n)
    end
    return x
end

function BR(game, A::Matrix{Float32}, B::Matrix{Float32}, V::Matrix{Float32}, V0::Matrix{Float32}, n::Int8)::Float32
    """
    Compute best reply of player n
    -------------
        # DEMAND    d := (A1 + A2) * exp(-s*x) / ( (A1 + A2) * exp(-s*x) + B )

        # PROFITS   p := d*(x - c) + d1*V1 + d2*V2 + d0*V0
            where   d1 := A1 * exp(- s*x) / ( A1*exp(-s*x) + B )

        # FOC       0 == e*B*(x - c) + (1-e)*(e1*V1 + e2*V2) - d*V0 - e*(e+B)/s
            where   e := (A1 + A2) * exp(- s*x)

        # BR        br = ( #1 + lambertw( exp(-#) * (A1 + A2) / B) ) / s
            where   # := 1 + s*(c + V0/B - (A1*V1 + A2*V2)/(A1 + A2) )
    """

    # Compute best reply
    arg = 1 + game.sigma*(game.c + V0[n]/B[n] - (A[n,1]*V[n,1] + A[n,2]*V[n,2])/(A[n,1] + A[n,2]) )
    br = ( arg + lambertw( exp(-arg) * (A[n,1] + A[n,2]) / B[n] )) / game.sigma
    return br
end

function compute_AB(game, x::Vector{Float32}, w::Matrix{Float32}, active::Vector{Int8}, Es::Vector{Float32})::Tuple{Matrix{Float32},Matrix{Float32},Matrix{Float32}}
    """Compute components for best reply function"""

    # Compute exp(price) of each product
    Ep = exp.(-game.sigma .* x) .* active

    # Compute p of each system
    p = prod([Es .* Ep; exp(2*game.sigma*game.v0)][game.outcomes_index], dims=1)

    # Compute A
    A = Es .* ((Es .* Ep)[game.BR_index[:,1:2]])

    # Compute B
    B = sum(p[game.nonoutcomes_index], dims=2) .* active

    # Compute V0
    V0 = sum(p[game.nonoutcomes_index] .* w'[game.nonoutcomes_index .+ [0;5;10;15]], dims=2)

    return A,B,V0
end

function compute_w_ownership(game, w, s, active_n)
    """Add extra components of w due to ownership"""
    o = game.ownership[s[5]+1,:]
    wown= Float32.(zeros(4,5)) + w
    for n in active_n
        if o[n]>0
            wown[n,:] += w[o[n],:]
        end
    end
    return wown
end

function compute_V_ownership(game, V, x, s, active_n)
    """Add extra components of w due to ownership"""
    o = game.ownership[s[5]+1,:]
    Vown = Float32.(zeros(4,2)) + V
    for n in active_n
        if o[n]>0
            Vown[n,1] += x[o[n]] - game.c
        end
    end
    return Vown
end

function update_p_BR(game, s::Vector{Int8}, w, p0::Vector{Float32})::Vector{Float32}
    """Update price in one state by BR iteration"""

    # Check active firms
    active = Int8.(s[1:4] .> 0)
    active_n = Int8.([n for n in 1:4 if active[n]>0])

    # Init x
    x = p0 .*  active[1:4]
    br = Float32.(zeros(4))

    # Compute exp(state) of each product
    Es = exp.(game.alpha .* s[1:4]) .* active

    # Modify w for ownership
    wown = compute_w_ownership(game, w, s, active_n)

    # Compute V
    V = wown'[game.outcomes_index[:,1:4]' .+ [0;5;10;15]]

    # Best until convergence
    dist = 1
    while dist>game.accuracy
        # Update BR
        for n in active_n
            A, B, V0 = compute_AB(game, x, wown, active, Es)
            Vown = compute_V_ownership(game, V, x, s, active_n)
            x[n] = BR(game, A, B, Vown, V0, n)
        end
        dist = max(abs.(x-br)...)
        br = copy(x)
    end
    return x
end

end
