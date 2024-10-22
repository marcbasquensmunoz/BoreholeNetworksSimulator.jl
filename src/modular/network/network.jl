

@with_kw struct BoreholeNetwork
    graph::SimpleGraph
end

struct Valve{T <: Number}
    split::Dict{Int, T}
end
act(v::Valve, i, m) = m * v.split[i]

@with_kw struct BoreholeOperation{T <: Number}
    valves::Dict{Int, Valve{T}} = Dict{Int, Valve{Float64}}()
    initial_mass_flow::T
end

# Maybe do something if a branch receives 0 mass flow
n_branches(n::BoreholeNetwork) = length(neighbors(n.graph, source(n.graph)))
BoreholeNetwork(n::Int) = BoreholeNetwork(graph=SimpleGraph(n+2))
source(network::BoreholeNetwork) = nv(network.graph) - 1
sink(network::BoreholeNetwork) = nv(network.graph)
connect!(network::BoreholeNetwork, i::Int, j::Int) = add_edge!(network.graph, i, j)
connect_to_source!(network::BoreholeNetwork, i) = add_edge!(network.graph, source(network), i)
connect_to_sink!(network::BoreholeNetwork, i) = add_edge!(network.graph, sink(network), i)
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
    visited = zeros(Bool, nv(network.graph))
    mass_flows .= zero(eltype(mass_flows))
    q = Queue{Tuple{Int, T}}()
    enqueue!(q, (source(network), operation.initial_mass_flow))
    mass_flows[source(network)] = operation.initial_mass_flow
    visited[source(network)] = true

    while !isempty(q)
        current_node, current_m = dequeue!(q)
        visited[current_node] = true
        if current_node in keys(operation.valves)
            valve = operation.valves[current_node]
            for n in neighbors(network.graph, current_node)
                if !visited[n]
                    new_m = act(valve, n, current_m)
                    mass_flows[n] += new_m
                    enqueue!(q, (n, new_m))
                end
            end
        elseif length(neighbors(network.graph, current_node)) < 2
            message = current_node == source(network) ? "source" : "$(current_node)"
            throw("Borehole $(message) splits into more than 1 branch; please provide a Valve to specify the flow going to each branch.")
        else
            for n in neighbors(network.graph, current_node)
                if !visited[n]
                    mass_flows[n] += current_m
                    enqueue!(q, (n, current_m))
                end
            end
        end
    end
end