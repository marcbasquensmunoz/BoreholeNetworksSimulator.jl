
struct SimpleOperator <: Operator
    mass_flows
end
SimpleOperator(;mass_flow, branches) = SimpleOperator(mass_flow .* ones(branches))

function operate(operator::SimpleOperator, i, options, Tin, Tout, Tb, q)
    network = options.configurations[1]
    BoreholeOperation(network, operator.mass_flows)
end
