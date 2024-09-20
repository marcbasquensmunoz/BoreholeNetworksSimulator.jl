using BoreholeNetworksSimulator
using BNSPlots

include("defs.jl")

Q_thr = 20000

Q_thr = 20000.
global toggle = 2
global single_branch = false
global hours_used = 0

mfs = zeros(2, Nt)

function operator(i, Tin, Tout, Tb, q, configurations, mass_flow_containers)
    mf = 0.5

    if Q_tot[i] > Q_thr
        mf1 = mf
        mf2 = mf
        global single_branch = false
        global hours_used = 0
    else
        if !single_branch || hours_used > 3
            global toggle = (toggle)%2 + 1
            global single_branch = true
            global hours_used = 0
        end
        mf1 = mf * (toggle%2)
        mf2 = mf * (toggle-1)%2
        global hours_used += 1
    end

    active_network = configurations[1]
    Nbr = n_branches(active_network)
    mass_flow_containers[1:6] .= mf1
    mass_flow_containers[7:10] .= mf2

    mfs[1, i] = mf1
    mfs[2, i] = mf2
    BoreholeOperation(active_network, @view mass_flow_containers[1:Nbr])
end

containers = @time initialize(options)
@time simulate!(operator=operator, options=options, containers=containers)

t_range = (5*8760-24*7):5*8760
toggle_plot = monitor(containers, [4, 7], t_range, options.t, color_pair = (colorant"darkgreen", colorant"red"), mf=mfs) 

CairoMakie.activate!()
save("examples/tekniska/plots/toggle.pdf", toggle_plot)