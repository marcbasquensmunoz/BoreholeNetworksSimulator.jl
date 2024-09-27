
"""
    SimpleOperator{T <: Number} <: Operator (
        mass_flows::Vector{T}
    )

Constant operation strategy, always returning `BoreholeOperation(options.configurations[1], operator.mass_flows)`.
"""
struct SimpleOperator{T <: Number} <: Operator
    mass_flows::Vector{T}
end
SimpleOperator(;mass_flow, branches) = SimpleOperator(mass_flow .* ones(eltype(mass_flow), branches))

function operate(operator::SimpleOperator, i, options, Tin, Tout, Tb, q)
    network = options.configurations[1]
    BoreholeOperation(network, operator.mass_flows)
end
