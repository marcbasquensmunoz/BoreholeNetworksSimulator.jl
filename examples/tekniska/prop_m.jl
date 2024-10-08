using BoreholeNetworksSimulator
using BNSPlots

include("defs.jl")

struct VariableMFOperator{T <: Number} <: Operator
    mass_flow_series::Vector{T}
    mass_flows::Vector{T}
end

function BoreholeNetworksSimulator.operate(operator::VariableMFOperator, step, options, Tin, Tout, Tb, q)
    operator.mass_flows .= operator.mass_flow_series[step]
    BoreholeOperation(options.configurations[1], operator.mass_flows)
end

operator = VariableMFOperator(0.5 .* Q_tot ./ Q_ref, zeros(n_branches(network)))

containers = @time initialize(options)
@time simulate!(operator=operator, options=options, containers=containers)

t_range = (5*8760-24*7):5*8760
prop_m_plot = monitor(containers, [4, 7], options.t, steps = t_range, color_pair = (colorant"darkgreen", colorant"red")) 

# save("examples/tekniska/plots/prop_m.png", prop_m_plot)
