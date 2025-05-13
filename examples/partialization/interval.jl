include("setup.jl")

struct IntervalicPartializationStrategy
    total_mass_flow
end

function BoreholeNetworksSimulator.operate(op::IntervalicPartializationStrategy, step, options, X)
    network = options.configurations[1]
    weights = floor(step / (24*30)) % 2 == 0 ? [0., 1.] : [1., 0.]
    source_valve = valve(Graphs.outneighbors(network.graph, source(network)), weights)
    BoreholeOperation(Dict(source(network) => source_valve), op.total_mass_flow, network)
end

operator = IntervalicPartializationStrategy(total_mass_flow)

reset!(options)
@time simulate!(operator=operator, options=options, containers=containers)
fig = monitor(containers, [1, 2], options.t, Δt = :month)
save("$(@__DIR__)/plots/interval_partialization.png", fig)
