using BoreholeNetworksSimulator
using BNSPlots
using Parameters

include("defs.jl")

@with_kw mutable struct ToggleOperator{T <: Number} <: Operator 
    Q_threshold::T
    toggle::Int = 2
    single_branch::Bool = false
    hours_used::Int = 0
    mass_flow::T
    mass_flow_containers::Vector{T}
end

function BoreholeNetworksSimulator.operate(operator::ToggleOperator, step, options, Tin, Tout, Tb, q)
    @unpack mass_flow, Q_threshold, single_branch, hours_used, toggle, mass_flow_containers = operator
    current_load = options.constraint.Q_tot[step]

    if current_load > Q_threshold
        mf1 = mass_flow
        mf2 = mass_flow
        single_branch = false
        hours_used = 0
    else
        if !single_branch || hours_used > 3
            toggle = (toggle)%2 + 1
            single_branch = true
            hours_used = 0
        end
        mf1 = mass_flow * (toggle%2)
        mf2 = mass_flow * (toggle-1)%2
        hours_used += 1
    end

    branch1 = [1, 6, 2, 7, 3, 8]
    branch2 = [4, 5, 9, 10]
    mass_flow_containers[branch1] .= mf1
    mass_flow_containers[branch2] .= mf2

    @pack! operator = single_branch, hours_used, toggle, mass_flow_containers
    BoreholeOperation(options.configurations[1], mass_flow_containers)
end

operator = ToggleOperator(Q_threshold = 20000., mass_flow = 0.5, mass_flow_containers = zeros(n_branches(network)))
containers = @time initialize(options)
@time simulate!(operator=operator, options=options, containers=containers)

t_range = (5*8760-24*7):5*8760
toggle_plot = monitor(containers, [4, 7], options.t, steps = t_range, color_pair = (colorant"darkgreen", colorant"red")) 

save("examples/tekniska/plots/toggle.png", toggle_plot)