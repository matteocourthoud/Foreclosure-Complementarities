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

# Column names
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

# Corresponding labels
const collabels = [
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

const marketnames = [
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

"""Get ylabels"""
function get_ylabels(policy)
    if policy=="nomergers"
        return figures.marketnames[1:3]
    elseif policy=="nobundling"
        return figures.marketnames[1:7]
    else
        return figures.marketnames
    end
end

"""Make all alluvial plots"""
function plot_alluvial(filenames)

    # Loop over policies and files
    for filename in filenames
        for policy in figures.policies

            # Import data
            df = DataFrame(CSV.File("output/transitions/$(figures.model)/$policy/$filename.csv"))
            xlabels = ["t=$t" for t in Int.([0; unique(df.t)])]
            ylabels = get_ylabels(policy)
            s0 = ones(length(ylabels)) ./ length(ylabels)
            Q = [Matrix(df[df.t.==1, 4:end]) for t in unique(df.t)]

            # Alluvial plot
            p = AlluvialPlot.alluvial(
                s0,
                Q,
                xlabels,
                ylabels,
                title = "State to State Transitions",
                dpi = 300,
                ytickfontsize = 6,
                xtickfontsize = 6,
            )

            # Save plot
            mkpath("output/alluvial/$(figures.model)/")
            png(p, "output/alluvial/$(figures.model)/$(policy)_$filename.png")
        end
    end

end




"""Make comparative statics plots"""
function get_stats(df, alphas, sigmas, col)
    v = zeros(length(sigmas), length(alphas))
    for (i, s) in enumerate(sigmas)
        for (j, a) in enumerate(alphas)
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
function make_heatmap(sumstats, alphas, sigmas, figtitle)
    # Plot Heatmap
    m = max(abs.(sumstats)...)
    p = contourf(
        sumstats,
        c = reverse(cgrad(:RdBu_11)),
        levels = 14,
        linewidth = 0,
        size = (550, 500),
        title = figtitle,
        xlabel = "Economies of Scale, α",
        xticks = (2:2:length(alphas), alphas[2:2:end]),
        ylabel = "Product Differentiation, σ",
        yticks = (2:2:length(sigmas), sigmas[2:2:end]),
        clims = (-m - 0.1, m + 0.1),
        dpi = 200,
        framestyle = :box,
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
            figtitle = figures.collabels[i]
            p = make_heatmap(sumstats, alphas, sigmas, figtitle)
            mkpath("output/compstats/$(figures.model)/$policy/")
            png(p, "output/compstats/$(figures.model)/$policy/$col.png")

            # Plot difference
            if policy != "baseline"
                df = DataFrame(
                    CSV.File("output/sumstats/$(figures.model)/baseline.csv"),
                )
                baseline = get_stats(df, alphas, sigmas, col)
                rel_sumstats = sumstats - baseline

                # Plot relative state
                figtitle = "Δ $(figures.collabels[i])"
                p = make_heatmap(rel_sumstats, alphas, sigmas, figtitle)
                png(p, "output/compstats/$(figures.model)/$policy/diff_$col.png")
            end

        end
    end
end

end
