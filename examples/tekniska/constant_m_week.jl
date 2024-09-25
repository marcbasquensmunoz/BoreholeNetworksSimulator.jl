using BoreholeNetworksSimulator
using BNSPlots
using CairoMakie

include("defs.jl")

mf = zeros(Nt)

function operator(i, Tin, Tout, Tb, q, configurations, mass_flow_containers)
    mf1 = .5
    mf2 = .5
    active_network = configurations[1]
    Nbr = n_branches(active_network)
    mass_flow_containers[1:6] .= mf1
    mass_flow_containers[7:10] .= mf2
    mf[i] = mf1
    BoreholeOperation(active_network, @view mass_flow_containers[1:Nbr])
end

containers = @time initialize(options)
@time simulate!(operator=operator, options=options, containers=containers)

t_range = (5*8760-24*7):5*8760
const_m_plot = monitor(containers, [4, 7], t_range, options.t, color_pair = (colorant"darkgreen", colorant"red"), mf=hcat(mf, mf)') 

save("examples/tekniska/plots/const_m.png", const_m_plot)


t_range = 1:Nt
const_m_plot_5_year = monitor(containers, [4, 7], t_range, options.t, color_pair = (colorant"darkgreen", colorant"red"), mf=hcat(mf, mf)') 

CairoMakie.activate!()
save("examples/tekniska/plots/const_m_5_years.png", const_m_plot_5_year)