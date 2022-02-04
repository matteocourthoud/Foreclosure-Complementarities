
module AlluvialPlot

using Plots
using ColorSchemes

"""
Make Alluvial plot

Inputs
------

s0: Vector
    Vector of initial states, dimension J (number of columns)

Q: Vector
    Vector of matrices, dimension I (number of rows)
    Each matrix should be of dimension I x J
"""
function alluvial(s0, Q, xlabels, ylabels; kwargs...)

    # Dimensions
    I = length(s0)
    J = length(Q)+1

    # State probabilities
    pr_s = zeros(I, J)
    pr_s[:, 1] = s0
    for j=2:J
        pr_s[:, j] = (pr_s[:, j-1]' * Q[j-1])'
    end

    # Bar coordinates
    ybars_bottom = zeros(size(pr_s));
    ybars_top = zeros(size(pr_s));
    for j=1:J
        ybars_bottom[:,j] = cumsum(pr_s[:,j]*0.9 .+ 0.1/(I-1)) .- 0.1/(I-1) - pr_s[:,j]*0.9;
        ybars_top[:,j] = cumsum(pr_s[:,j]*0.9 .+ 0.1/(I-1)) .- 0.1/(I-1);
    end

    # Init plot
    colors = cgrad(:lajolla, I+2, categorical = true)[2:end-1]
    ypos = (ybars_top[:,1] + ybars_bottom[:,1])/2
    p = plot(xticks=(1:I+1, xlabels),
             yticks=(ypos, ylabels),
             axis=false,
             legend=false,
             grid=false;
             kwargs...)
    # Setup
    L = 100
    c = (1 .- cos.(range(0, pi, length=L)))./2;
    w = size(ybars_bottom,2)/40;

    # Plot shapes
    for i=1:I
        for j=1:J
            x_corners = [j-w, j+w, j+w, j-w];
            y_corners = [ybars_bottom[i,j], ybars_bottom[i,j], ybars_top[i,j], ybars_top[i,j]];
            plot!(p, Shape(x_corners, y_corners), color=colors[i], linealpha=0)
        end
    end

    # Add curve
    for j=1:J-1
        bottom_rights = ybars_bottom[:,j+1];
        for i=1:I

            # Get corners
            bottom_lefts = (cumsum(Q[j][i,:]) - Q[j][i,:]) * pr_s[i,j]*0.9 .+ ybars_bottom[i,j]
            top_lefts = cumsum(Q[j][i,:]) * pr_s[i,j]*0.9 .+ ybars_bottom[i,j]
            top_rights = bottom_rights + top_lefts - bottom_lefts

            # Get xs
            x1 = j+w
            x2 = j+1-w

            # Plot shapes
            for k=1:I
                y = [bottom_lefts[k] .+ c * (bottom_rights[k] - bottom_lefts[k]); top_rights[k] .+ c * (top_lefts[k] - top_rights[k])]
                x = [range(x1, x2, length=L); range(x2, x1, length=L)]
                plot!(p, Shape(x, y), color=colors[i], alpha=0.3, linealpha=0)
            end

            bottom_rights = top_rights
        end
    end


    return p
end

end # module
