using BoreholeNetworksSimulator
using BNSPlots

include("defs.jl")

operator = SimpleOperator(mass_flow=0.5, branches=n_branches(network))
containers = @time initialize(options)
@time simulate!(operator=operator, options=options, containers=containers)

t_range = (5*8760-24*7):5*8760
const_m_plot = monitor(containers, [4, 7], t_range, options.t, color_pair = (colorant"darkgreen", colorant"red")) 
save("examples/tekniska/plots/const_m.png", const_m_plot)


t_range = 1:Nt
const_m_plot_5_year = monitor(containers, [4, 7], t_range, options.t, color_pair = (colorant"darkgreen", colorant"red"), mf=hcat(mf, mf)') 
save("examples/tekniska/plots/const_m_5_years.png", const_m_plot_5_year)