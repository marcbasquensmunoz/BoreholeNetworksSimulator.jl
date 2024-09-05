using WGLMakie

function monitor_branch(containers::SimulationContainers, branch, color_pair, t; Δt = 8760/12.)
    scene = Figure(size = (600, 450))

    axis_T = scene[1, 1] = Axis(scene, ylabel = "T [°C]") #, xlabel = "time [months]")
    axis_Q = scene[2, 1] = Axis(scene, ylabel = "q [W/m]", xlabel = "time [years]")

    color_range = make_color_range(color_pair, length(branch)) 

    Tfin = BoreholeNetworksSimulator.Tfin(containers, branch)   
    q = BoreholeNetworksSimulator.q(containers, branch)

    @show Tfin

    for (borehole, color) in enumerate(color_range)
        lines!(axis_T, collect(t ./ 12Δt), Tfin[borehole, :], color = color, linewidth = 2.)
        lines!(axis_Q, collect(t ./ 12Δt), q[borehole, :], color = color, linewidth = 2.)
    end
    #scatter!(axis_T, collect(options.t ./ 12Δt), Tfos[:,nbranch1],  color = :black, markersize = 8.)
    #scatter!(axis_T, collect(options.t ./ 12Δt), Tfo,  color = :grey, markersize = 12.)

    hidexdecorations!(axis_Q, grid = false)

    scene
end