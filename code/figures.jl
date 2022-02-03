"""Figures and plots"""

module figures

using Plots
using CSV
using DataFrames

include("AlluvialPlot/src/AlluvialPlot.jl")

# Model
const model = "privacy"

# Policies
const policies = [
    "baseline",
    "nolearning",
    "nomergers",
    "nobundling",
    "datasharing",
    "limitedmergers",
    "limitedbundling",
    "nopredexitpricing",
    "nopredexitbundling",
    "nopredentrypricing",
    "nopredentrybundling",
]

const colnames = [
    "margin",
    "bcost",
    "entry",
    "exit",
    "merger",
    "bundling",
    "profits",
    "surplus",
    "welfare",
    "mpoly",
    "mpolyA",
    "mpolyB",
    "profitsLR",
    "surplusLR",
    "welfareLR",
]

const varnames = [
    "Price - Cost (short run)",
    "Below Cost Pricing (short run)",
    "Entry Probability (short run)",
    "Exit Probability (short run)",
    "Merger Probability (short run)",
    "Bundling Probability (short run)",
    "Total Profits (NPV)",
    "Consumer Surplus (NPV)",
    "Total Welfare (NPV)",
    "Monopoly Probability (long run)",
    "Monopoly Probability A (long run)",
    "Monopoly Probability B (long run)",
    "Total Profits (long run)",
    "Consumer Surplus (long run)",
    "Total Welfare (long run)",
]

const matketnames = [
    "Monopoly in A&B",
    "M Monopoly in A&B",
    "B Monopoly in A&B",
    "mixed M/Dpoly",
    "M mixed M/Dpoly",
    "B mixed M/Dpoly",
    "Duopoly in A&B",
    "M/2 Duopoly in A&B",
    "M Duopoly in A&B",
    "B/2 Duopoly in A&B",
    "B Duopoly in A&B",
]

const filenames = ["a80g0s80", "a30g0s70", "a70g0s30", "a70g0s70"]

"""Make alluvial plots"""
function plot_alluvial()
    # Input
    ylabels = ("state 1", "state 2", "state 3")
    xlabels = ("time 1", "time 2", "time 3", "time 4")
    s0 = [0.1, 0.2, 0.7]
    Q = [i ./ sum(i, dims = 2) for i in [rand(3, 3) for j = 1:3]]

    # Alluvial plot
    p = AlluvialPlot.alluvial(
        s0,
        Q,
        xlabels = xlabels,
        ylabels = ylabels,
        title = "Test 1",
    )

end

"""Make comparative statics plots"""
function get_stats(df, alphas, sigmas, col)
    v = zeros(length(alphas), length(sigmas))
    for (i, a) in enumerate(alphas)
        for (j, s) in enumerate(sigmas)
            v[i, j] = df[
                (df.alpha.==a).&
                (df.beta.==0.95).&
                (df.gamma.==0.0).&
                (df.p0.==1).&
                (df.market.==400).&
                (df.sigma.==s),
                col,
            ][1]
        end
    end
    return v
end

"""Plot heatmap"""
function make_heatmap(sumstats, alphas, sigmas, i)
    # Plot Heatmap
    m = max(abs.(sumstats)...)
    p = contourf(
        sumstats,
        c = :balance,
        levels = 10,
        linewidth = 0,
        size = (550, 500),
        title = figures.varnames[i],
        xlabel = "Product Differentiation, σ",
        xticks = (2:2:length(sigmas), sigmas[2:2:end]),
        ylabel = "Economies of Scale, α",
        yticks = (2:2:length(alphas), alphas[2:2:end]),
        clims = (-m - 0.1, m + 0.1),
        dpi = 300,
    )
end

"""Make comparative statics plots"""
function plot_compsats()

    # Loop over policies and statistics
    for policy in figures.policies
        for i = 1:length(figures.colnames)
            col = figures.colnames[i]
            print("\nPlotting $policy $col ($i)")

            # Import data
            df = DataFrame(CSV.File("output/sumstats/$(figures.model)/$policy.csv"))
            alphas = sort(unique(df.alpha))
            sigmas = sort(unique(df.sigma))
            sumstats = get_stats(df, alphas, sigmas, col)

            # Plot baseline
            p = make_heatmap(sumstats, alphas, sigmas, i)
            mkpath("output/compstats/$(figures.model)/$policy/")
            png(p, "output/compstats/$(figures.model)/$policy/$col.png")

            # Plot difference
            if policy != "baseline"
                df = DataFrame(
                    CSV.File("output/sumstats/$(figures.model)/baseline.csv"),
                )
                baseline = get_stats(df, alphas, sigmas, col)
                rel_sumstats = sumstats - baseline
                p = make_heatmap(rel_sumstats, alphas, sigmas, i)
                png(p, "output/compstats/$(figures.model)/$policy/diff_$col.png")
            end

        end
    end
end

end
