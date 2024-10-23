using Graphs

make_color_range(color_pair, n) = n == 1 ? [color_pair[1]] : range(color_pair[1], stop=color_pair[2], length=n)

"""
    plot_borefield(network, positions; distinguished_boreholes = [])

Makes a plot of the borefield, showing the boreholes numbered and their connections.
# Arguments
- `network`: Network specifying the connections between boreholes.
- `positions`: The positions of each borehole.

# Optional arguments
- `distinguished_boreholes`: Vector of `Tuple{Int, Color}`. If specified, the boreholes corresponding to the given values will be highlighted with each of the colors provided.
"""
function plot_borefield(network, positions; distinguished_boreholes = [])
    scene = Figure()
    axis = scene[1, 1] = Axis(scene, ylabel = "y [m]", xlabel = "x[m]", aspect = DataAspect())
    min_x = minimum(map(x->x[1], positions))
    max_x = maximum(map(x->x[1], positions))
    min_y = minimum(map(x->x[2], positions))
    max_y = maximum(map(x->x[2], positions))

    margin = max(max_x-min_x, max_y-min_y) * 0.1

    xlims!(axis, min_x-margin, max_x+margin)
    ylims!(axis, min_y-margin, max_y+margin)

    N = nv(network.graph)
    borefield = SimpleGraph(network.graph)
    rem_vertex!(borefield, N)
    rem_vertex!(borefield, N - 1)

    Nb = nv(borefield)
    borehole_colors = [:black for i in 1:Nb]
    borehole_sizes = 12 * ones(Int, Nb)
    for (bh, color) in distinguished_boreholes
        borehole_colors[bh] = color
        borehole_sizes[bh] = 16
    end
    graphplot!(axis, borefield, layout=positions, nlabels=["$i" for i in 1:Nb], node_color = borehole_colors, node_size = borehole_sizes, nlabels_distance = 5.)

    hidespines!(axis)

    scene   
end