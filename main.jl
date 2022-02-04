#=
Foreclosure Complementarities
By: Matteo Courthoud
Date: 08/06/2021
=#

include("code/solve.jl")
include("code/figures.jl")

# Parameters
model = "privacy"
alphas = collect(0.0:0.1:1.0)
gammas = collect(0.0:0.1:0.0)
sigmas = collect(1.0:1.0:10.0)
policies = ["baseline",
            "nolearning",
            "nomergers",
            "nobundling",
            "datasharing",
            "limitedmergers",
            "limitedbundling",
            "nopredexitpricing",
            "nopredexitbundling",
            "nopredentrypricing",
            "nopredentrybundling"]

# Replicate
solve.replicate(model, alphas, gammas, sigmas, policies)

# Make graphs
figures.plot_compsats()
figures.plot_alluvial(("a80g0s80",))
