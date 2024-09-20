
make_color_range(color_pair, n) = n == 1 ? [color_pair[1]] : range(color_pair[1], stop=color_pair[2], length=n)

"""
    plot_borefield(network, positions; distinguished_branches = [], colors = [])

Makes a plot of the borefield, showing the boreholes numbered and their connections.
# Arguments
- `network`: Network specifying the connections between boreholes.
- `positions`: The positions of each borehole.

# Optional arguments
- `distinguished_branches`: Vector of `Int`. If specified, the branches corresponding to the given values will be highlighted with the colors of `colors`.
- `colors`: Vector of `Tuples` of `Color`. If specified, uses each pair of colors to define a color gradient to color each of the branches of `distinguished_branches`.
"""
function plot_borefield(network, positions; distinguished_branches = [], colors = [])

    scene = Figure()
    axis = scene[1, 1] = Axis(scene, ylabel = "y [m]", xlabel = "x[m]", aspect = DataAspect())
    min_x = minimum(map(x->x[1], positions))
    max_x = maximum(map(x->x[1], positions))
    min_y = minimum(map(x->x[2], positions))
    max_y = maximum(map(x->x[2], positions))
    
    margin = max(max_x-min_x, max_y-min_y) * 0.1

    xlims!(axis, min_x-margin, max_x+margin)
    ylims!(axis, min_y-margin, max_y+margin)

    scatter!(axis,  Makie.Point2.(positions), markersize = 15.)

    # Draw boreholes
    for (i, n_branch) in enumerate(distinguished_branches)
        branch = network.branches[n_branch]
        branch_colors = make_color_range(colors[i], length(branch))
        for (borehole, color) in zip(branch, branch_colors)
            point =  [ Makie.Point(positions[borehole]) ]
            scatter!(axis, point, color = color, markersize = 24)
        end
    end

    # Write labels with borehole numbers
    for i in eachindex(1:length(positions))
        text!(axis, "$i", fontsize = 9, position  = (positions[i] .+  (0.7 , -0.7)) , color = :blue)
    end

    # Draw lines representing connections
    for branch in network.branches
        for k in 2:length(branch)
            s = Makie.Point.([ positions[branch][k-1], positions[branch][k] ] )
            linesegments!(axis, s, color = :red, linewidth = 1.5)
        end
    end

    hidespines!(axis)

    scene   
end