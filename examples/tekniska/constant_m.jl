using BoreholeNetworksSimulator
using BNSPlots

include("defs.jl")

operator = ConstantOperator(network, mass_flows=0.5 * ones(10))
containers = initialize(options)
simulate!(operator=operator, options=options, containers=containers)

t_range = (5*8760-24*7):5*8760
const_m_plot = monitor(containers, [4, 7], options.t, steps = t_range, colors = [colorant"darkgreen", colorant"red"]) 
# save("examples/tekniska/plots/const_m.png", const_m_plot)

const_m_plot_5_year = monitor(containers, [4, 7], options.t, colors = [colorant"darkgreen", colorant"red"]) 
# save("examples/tekniska/plots/const_m_5_years.png", const_m_plot_5_year)