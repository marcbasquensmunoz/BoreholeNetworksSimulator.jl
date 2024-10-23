

@with_kw struct BoreholeNetwork
    graph::SimpleDiGraph = SimpleDiGraph()
end

struct Valve{T <: Number}
    split::Dict{Int, T}
end
act(v::Valve, i, m) = m * v.split[i]

struct BoreholeOperation{T <: Number}
    valves::Dict{Int, Valve{T}}
    source_mass_flow::T
    network::BoreholeNetwork
end

BoreholeOperation(::Nothing) = BoreholeOperation(Dict{Int, Valve{Float64}}(), 0., BoreholeNetwork())
BoreholeOperation(;network::BoreholeNetwork, source_mass_flow::T) where {T <: Number} = BoreholeOperation(Dict{Int, Valve{T}}(), source_mass_flow, network)
function BoreholeOperation(;network::BoreholeNetwork, mass_flows::Vector{T}) where {T <: Number}
    source_valve = absolute_valve(outneighbors(network.graph, source(network)), mass_flows)
    BoreholeOperation(Dict(source(network) => source_valve), sum(mass_flows), network)
end
# Maybe do something if a branch receives 0 mass flow
n_boreholes(n::BoreholeNetwork) = nv(n.graph) - 2
n_branches(n::BoreholeNetwork) = length(neighbors(n.graph, source(n)))
initialize_mass_flows(network::BoreholeNetwork) = zeros(nv(network.graph))
first_bhs_in_branch(network::BoreholeNetwork) = outneighbors(network.graph, source(network))
BoreholeNetwork(n::Int) = BoreholeNetwork(graph=SimpleDiGraph(n+2))
source(network::BoreholeNetwork) = nv(network.graph) - 1
sink(network::BoreholeNetwork) = nv(network.graph)
connect!(network::BoreholeNetwork, i::Int, j::Int) = add_edge!(network.graph, i, j)
connect_to_source!(network::BoreholeNetwork, i::Int) = add_edge!(network.graph, source(network), i)
connect_to_source!(network::BoreholeNetwork, I::Vector{Int}) = foreach(i -> add_edge!(network.graph, source(network), i), I)
connect_to_sink!(network::BoreholeNetwork, i::Int) = add_edge!(network.graph, i, sink(network))
connect_to_sink!(network::BoreholeNetwork, I::Vector{Int}) = foreach(i -> add_edge!(network.graph, i, sink(network)), I)
function connect_in_series!(network::BoreholeNetwork, I::Vector{Int}) 
    for i in 2:length(I)
        connect!(network, I[i-1], I[i])
    end
end
function connect_in_parallel!(network::BoreholeNetwork, j::Int, I::Vector{Int}) 
    for i in I
        connect!(network, j, i)
    end
end
function equal_valve(nodes::Vector{Int}) 
    valve = Valve(Dict{Int, Float64}())
    for n in nodes
        valve.split[n] = 1/length(nodes)
    end
    return valve
end
function valve(nodes::Vector{Int}, weights::Vector{T}) where {T <: Number}
    valve = Valve(Dict{Int, T}())
    for i in 1:length(nodes)
        valve.split[nodes[i]] = weights[i] 
    end
    return valve
end

function absolute_valve(nodes::Vector{Int}, mass_flows::Vector{T}) where {T <: Number}
    valve = Valve(Dict{Int, T}())
    total_mass_flow = sum(mass_flows)
    for (i, n) in enumerate(nodes)
        valve.split[n] = mass_flows[i] / total_mass_flow
    end
    return valve
end

function all_parallel_network(n) 
    network = BoreholeNetwork(n)
    for i in 1:n
        connect_to_source!(network, i)
        connect_to_sink!(network, i)
    end
    return network
end

function all_series_network(n) 
    network = BoreholeNetwork(n)
    connect_to_source!(network, 1)
    for i in 2:n
        connect!(network, i-1, i)
    end
    connect_to_sink!(network, n)
    return network
end

function compute_mass_flows!(mass_flows, network::BoreholeNetwork, operation::BoreholeOperation{T}) where {T <: Number}
    mass_flows .= zero(eltype(mass_flows))
    q = Queue{Tuple{Int, T}}()
    enqueue!(q, (source(network), operation.source_mass_flow))
    mass_flows[source(network)] = operation.source_mass_flow

    while !isempty(q)
        current_node, current_m = dequeue!(q)
        if current_node in keys(operation.valves)
            valve = operation.valves[current_node]
            for n in outneighbors(network.graph, current_node)
                new_m = act(valve, n, current_m)
                mass_flows[n] += new_m
                enqueue!(q, (n, new_m))
            end
        elseif length(outneighbors(network.graph, current_node)) > 1
            message = current_node == source(network) ? "source" : "$(current_node)"
            throw("Borehole $(message) splits into more than 1 branch; please provide a Valve to specify the flow going to each branch.")
        else
            for n in outneighbors(network.graph, current_node)
                mass_flows[n] += current_m
                enqueue!(q, (n, current_m))
            end
        end
    end
end