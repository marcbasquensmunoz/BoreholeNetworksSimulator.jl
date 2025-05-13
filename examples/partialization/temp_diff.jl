using Parameters
using BNSPlots: get_Tb_range
include("setup.jl")

@with_kw mutable struct TemperatureDifferencePartializationStrategy{T <: Number} <: Operator 
    current_bh::Int = 1
    Tb_diff::T
    mass_flow::T
    mass_flow_containers::Vector{T}
end

function BoreholeNetworksSimulator.operate(operator::TemperatureDifferencePartializationStrategy, step, options, X)
    @unpack mass_flow, Tb_diff, mass_flow_containers, current_bh = operator

    if step != 1    
        @views Tb = X[get_Tb_range(2), :]
        if abs(Tb[1, step-1] - Tb[2, step-1]) > Tb_diff
            current_bh = current_bh%2 + 1
        end
    end

    mf1 = 0.
    mf2 = 0.
    if current_bh == 1
        mf1 = mass_flow
    elseif current_bh == 2
        mf2 = mass_flow
    end
    mass_flow_containers .= [mf1, mf2]
    @pack! operator = current_bh
    BoreholeOperation(network = options.configurations[1], mass_flows = mass_flow_containers)
end

temp_diff_operator = TemperatureDifferencePartializationStrategy(
    Tb_diff = .2,
    mass_flow = total_mass_flow, 
    mass_flow_containers = zeros(n_branches(network))
)

reset!(options)
@time simulate!(operator=temp_diff_operator, options=options, containers=containers)
fig = monitor(containers, [1, 2], options.t, Δt = :month)
save("$(@__DIR__)/plots/temp_diff_partialization.png", fig)