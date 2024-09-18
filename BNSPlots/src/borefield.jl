
make_color_range(color_pair, n) = range(color_pair[1], stop=color_pair[2], length=n)

function plot_borefield(network, positions; distinguished_branches = [], colors = [], size = (400, 400))

    scene = Figure(size=size)
    axis = scene[1, 1] = Axis(scene, ylabel = "y [m]", xlabel = "x[m]", aspect = DataAspect())
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
        text!(axis, "$i", fontsize = 9, position  = (positions[i] .+  (0.8 , -0.7)) , color = :blue)
    end

    # Draw lines representing connections
    for branch in network.branches
        for k in 2:length(branch)
            s = Makie.Point.([ positions[branch][k-1], positions[branch][k] ] )
            linesegments!(axis, s, color = :red, linewidth = 1.5)
        end
    end

    scene   
end