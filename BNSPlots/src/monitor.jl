
function anyin(elements, v)
    for element in elements
        if element in v
            return true
        end
    end
    return false
end

"""
    monitor(containers, branch, t; display = [:Tfin, :Tfout, :Tb, :q], Δt = :year, color_pair = (colorant"navajowhite2", colorant"darkgreen")

Creates a plot of the result of the simulation.
# Arguments
- `containers`: The containers (`SimulationContainers`) containing the result of the simulation through `simulate!`.
- `branch`: A vector containing the IDs of the boreholes whose data will be displayed. 
- `t`: The times at which the data corresponds. It should normally be `options.t`.

# Optional arguments
- `display`: A vector describing which plots that will be generated. If `:Tfin`, `:Tfout` or `:Tb` are specified, a temperature plot will be created showing the inlet fluid temperature, 
    the outlet fluid temperature, and the borehole wall temperature, respectively. If `:q` is specified, a separate power plot will be created shwoing the heat extracted per meter.
- `Δt`: The scale of the x-axis in the plot. Possible options: :year, :month, :hour.
- `color_pair`: A pair of colors used as extrema to generate a range of colors for each borehole. 
"""
function monitor(containers, branch, steps, t; display = [:Tfin, :Tfout, :Tb, :q, :mf], Δt = :year, color_pair = (colorant"navajowhite2", colorant"darkgreen"), mf)
    if isempty(display)
        return
    end
    if !(Δt in [:year, :month, :hour])
        return throw("Plot time step must be one of :year, :month, :hour")
    end

    scene = Figure(size=(800, 500))
    grid = scene[1, 1] = GridLayout()
    axes = []

    color_range = make_color_range(color_pair, length(branch)) 
    
    if anyin([:Tfin, :Tfout, :Tb], display)
        axis_T = Axis(grid[length(axes)+1, 1], ylabel = L"T \, \left[ °C \right]")
        push!(axes, axis_T)
    end
    if :q in display
        axis_Q = Axis(grid[length(axes)+1, 1], ylabel = L"q \, \left[ \frac{W}{m} \right]")
        push!(axes, axis_Q)
    end
    if :mf in display
        axis_m = Axis(grid[length(axes)+1, 1], ylabel = L"\dot{m} \, \left[ \frac{kg}{s} \right]")
        push!(axes, axis_m)
    end

    Tfin = get_Tfin(containers, branch)  
    Tfout = get_Tfout(containers, branch)  
    Tb = get_Tb(containers, branch)    
    q = get_q(containers, branch)

    secs_in_year = 8760*3600
    conversion = Dict(:year => 1, :month => 12, :hour => 8760)
    time = collect(t[steps] ./ secs_in_year * conversion[Δt])

    for (borehole, color) in enumerate(color_range)
        if :Tfin in display
            lines!(axis_T, time, Tfin[borehole, steps], color = color, linewidth = 2.)
        end
        if :Tfout in display
            lines!(axis_T, time, Tfout[borehole, steps], color = color, linewidth = 1.2, linestyle=:dash)
        end
        if :Tb in display
            scatter!(axis_T, time, Tb[borehole, steps], color = color, markersize = 5.)
        end
        if :q in display
            lines!(axis_Q, time, q[borehole, steps], color = color, linewidth = 2.)
        end
        if :mf in display
            stairs!(axis_m, time, mf[borehole, steps], color = color, linewidth = 2.)
        end
    end


    group_color = [PolyElement(color = color, strokecolor = :transparent) for color in color_range]
    group_marker = [LineElement(color = :black), LineElement(color = :black, linestyle = :dash), MarkerElement(marker = :circle, color = :black, strokecolor = :transparent, markersize = 5.)]
    legend = Legend(scene, [group_color, group_marker], [string.(branch), ["Tfin", "Tfout", "Tb"]], ["Boreholes", "Temperatures"], tellheight = true, tellwidth = false)
    legend.titleposition = :top
    legend.orientation = :horizontal
    legend.nbanks = 1
    scene[2,1] = legend
    for i in 1:length(axes)-1
        hidexdecorations!(axes[i], grid = false)
        linkxaxes!(axes[i], axes[i+1])
    end
    axes[end].xlabel = "time [$(Δt)s]"
    rowgap!(grid, 10)

    scene
end