
function anyin(elements, v)
    for element in elements
        if element in v
            return true
        end
    end
    return false
end

"""
    monitor(containers, branch, t; steps = 1:length(t), display = [:Tfin, :Tfout, :Tb, :q], Δt = :year, color_pair = (colorant"navajowhite2", colorant"darkgreen")

Creates a plot of the result of the simulation.
# Arguments
- `containers`: The containers (`SimulationContainers`) containing the result of the simulation through `simulate!`.
- `boreholes`: A vector containing the IDs of the boreholes whose data will be displayed. 
- `t`: The times at which the data corresponds. It should normally be `options.t`.

# Optional arguments
- `steps`: Index of the steps to display in the plot.
- `display`: A vector describing which plots that will be generated. If `:Tfin`, `:Tfout` or `:Tb` are specified, a temperature plot will be created showing the inlet fluid temperature, 
    the outlet fluid temperature, and the borehole wall temperature, respectively. If `:q` is specified, a separate power plot will be created shwoing the heat extracted per meter.
- `Δt`: The scale of the x-axis in the plot. Possible options: `:year`, `:month`, `:hour`.
- `colors`: A list of colors used for each borehole. If not specified, the colors used will be between colorant"navajowhite2" and colorant"darkgreen".
"""
function monitor(containers, boreholes, t; steps = 1:length(t), display = [:Tfin, :Tfout, :Tb, :q, :mf], Δt = :year, colors = [])
    if isempty(display)
        return
    end
    if !(Δt in [:year, :month, :hour])
        return throw("Plot time step must be one of :year, :month, :hour")
    end

    scene = Figure(size=(800, 500))
    grid = scene[1, 1] = GridLayout()
    axes = []

    if !isempty(colors)
        color_range = colors
    else
        color_range = make_color_range((colorant"navajowhite2", colorant"darkgreen"), length(boreholes)) 
    end

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

    Tfin = get_Tfin(containers, boreholes)  
    Tfout = get_Tfout(containers, boreholes)  
    Tb = get_Tb(containers, boreholes)    
    q = get_q(containers, boreholes)
    mf = get_mf(containers, boreholes)

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
    group_marker = []
    legend_markers = []
    if :Tfin in display
        push!(legend_markers, L"T_{f, \text{in}}")
        push!(group_marker, LineElement(color = :black, linestyle = :solid))
    end
    if :Tfout in display
        push!(legend_markers, L"T_{f, \text{out}}")
        push!(group_marker, LineElement(color = :black, linestyle = :dash))
    end
    if :Tb in display
        push!(legend_markers, L"T_b")
        push!(group_marker, MarkerElement(marker = :circle, color = :black, strokecolor = :transparent, markersize = 5.))
    end
    legend = Legend(scene, 
        [group_color, group_marker], [string.(boreholes), legend_markers], ["Boreholes", "Temperatures"], 
        tellheight = true, 
        tellwidth = false
    )
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