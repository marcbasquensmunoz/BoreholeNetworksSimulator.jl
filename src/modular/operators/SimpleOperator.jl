
"""
    ConstantOperator{T <: Number} <: Operator (
        valves::Dict{Int, Valve{T}} = Dict{Int, Valve{Float64}}()
        mass_flow::T
    )

Constant operation strategy. Its implementation of `operate` is:
`operate(op::ConstantOperator, step, options, X) = BoreholeOperation(op.valves, op.mass_flow, options.configurations[1])`

# Convenience constructors

    ConstantOperator(network::BoreholeNetwork; mass_flows) 

Builds the initial valve and total source mass flow needed to have mass flow equal to `mass_flows[i]` in branch `i`, according to the network `network`.
"""
@with_kw struct ConstantOperator{T <: Number} <: Operator
    valves::Dict{Int, Valve{T}} = Dict{Int, Valve{Float64}}()
    mass_flow::T
end

function ConstantOperator(network::BoreholeNetwork; mass_flows) 
    source_valve = absolute_valve(outneighbors(network.graph, source(network)), mass_flows)
    ConstantOperator(valves = Dict(source(network) => source_valve), mass_flow = sum(mass_flows))
end
operate(op::ConstantOperator, step, options, X) = BoreholeOperation(op.valves, op.mass_flow, options.configurations[1])
