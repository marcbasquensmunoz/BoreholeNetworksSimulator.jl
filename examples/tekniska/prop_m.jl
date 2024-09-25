using BoreholeNetworksSimulator
using BNSPlots

include("defs.jl")

mf = zeros(Nt)

function operator(i, Tin, Tout, Tb, q, configurations, mass_flow_containers)
    m = 0.5 * Q_tot[i]/Q_ref
    mf1 = m
    mf2 = m
    active_network = configurations[1]
    Nbr = n_branches(active_network)
    mass_flow_containers[1:6] .= mf1
    mass_flow_containers[7:10] .= mf2
    mf[i] = m
    BoreholeOperation(active_network, @view mass_flow_containers[1:Nbr])
end

@with_kw struct VariableMFOperator <: Operator
    mass_flows
end

function operate(operator::VariableMFOperator, step, options, Tin, Tout, Tb, q)

    BoreholeOperation(active_network, @view mass_flow_containers[1:Nbr])
end

operator = VariableMFOperator(mass_flows = 0.5 * Q_tot[i]/Q_ref)

containers = @time initialize(options)
@time simulate!(operator=operator, options=options, containers=containers)

t_range = (5*8760-24*7):5*8760
prop_m_plot = monitor(containers, [4, 7], t_range, options.t, color_pair = (colorant"darkgreen", colorant"red"), mf=hcat(mf, mf)') 

save("examples/tekniska/plots/prop_m.png", prop_m_plot)