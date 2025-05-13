
"""
    BoreholeNetwork(
        graph::SimpleDiGraph = SimpleDiGraph()
    )

Represents the network of boreholes with a directed graph `graph`, whose edges represent the flow of fluid.
The `graph` is expected to have number of vertices `Nb+2`, where `Nb` is the number of boreholes in the borefield.
Vertexs number `Nb+1` and `Nb+2` represent the mass flow source and sink, respectively. 
All branches are expected to be connected to the source at the beggining and to the sink at the end.

# Convenience constructors
 - `BoreholeNetwork(n::Int)`
Create a network of `n` boreholes (and the source and sink) without any connections.

 - `all_parallel_network(n::Int)`
Create a network of `n` boreholes connected in parallel.

 - `all_series_network(n::Int)`
Create a network of `n` boreholes connected in series.
"""
@with_kw struct BoreholeNetwork
    graph::SimpleDiGraph{Int} = SimpleDiGraph{Int}()
end

BoreholeNetwork(n::Int) = BoreholeNetwork(graph=SimpleDiGraph(n+2))
n_boreholes(n::BoreholeNetwork) = nv(n.graph) - 2
n_branches(n::BoreholeNetwork) = length(neighbors(n.graph, source(n)))
initialize_mass_flows(network::BoreholeNetwork) = zeros(nv(network.graph))
"""
    first_bhs_in_branch(network::BoreholeNetwork)

Return the first borehole in each branch of the borefield.
"""
first_bhs_in_branch(network::BoreholeNetwork) = outneighbors(network.graph, source(network))

"""
    source(network::BoreholeNetwork) 

Return the vertex index of the source.
"""
source(network::BoreholeNetwork) = nv(network.graph) - 1
"""
    sink(network::BoreholeNetwork) 

Return the vertex index of the sink.
"""
sink(network::BoreholeNetwork) = nv(network.graph)

"""
    connect!(network::BoreholeNetwork, i::Int, j::Int)

Create a flow from borehole `i` to borehole `j`.
"""
connect!(network::BoreholeNetwork, i::Int, j::Int) = add_edge!(network.graph, i, j)
"""
    connect_to_source!(network::BoreholeNetwork, i::Int)

Create a connection injecting mass flow to borehole `i`.
"""
connect_to_source!(network::BoreholeNetwork, i::Int) = add_edge!(network.graph, source(network), i)
"""
    connect_to_source!(network::BoreholeNetwork, I::Vector{Int})

Vector version of connect_to_source!.
"""
connect_to_source!(network::BoreholeNetwork, I::Vector{Int}) = foreach(i -> add_edge!(network.graph, source(network), i), I)
"""
    connect_to_sink!(network::BoreholeNetwork, i::Int)

Create a connection extracting mass flow from borehole `i`.
"""
connect_to_sink!(network::BoreholeNetwork, i::Int) = add_edge!(network.graph, i, sink(network))
"""
    connect_to_sink!(network::BoreholeNetwork, I::Vector{Int})

Vector version of connect_to_sink!.
"""
connect_to_sink!(network::BoreholeNetwork, I::Vector{Int}) = foreach(i -> add_edge!(network.graph, i, sink(network)), I)
"""
    connect_in_series!(network::BoreholeNetwork, I::Vector{Int}) 

Connect all the boreholes in `I`, in order.
"""
function connect_in_series!(network::BoreholeNetwork, I::Vector{Int}) 
    for i in 2:length(I)
        connect!(network, I[i-1], I[i])
    end
end
"""
    connect_in_parallel!(network::BoreholeNetwork, j::Int, I::Vector{Int}) 

Make a connection from borehole `j` to each of the boreholes in `I`.
"""
function connect_in_parallel!(network::BoreholeNetwork, j::Int, I::Vector{Int}) 
    for i in I
        connect!(network, j, i)
    end
end
"""
    boreholes_in_branch(n::BoreholeNetwork; first_bh::Int)

Return all the boreholes in the branch that starts by borehole `first_bh`.
"""
function boreholes_in_branch(n::BoreholeNetwork; first_bh::Int)
    boreholes = neighborhood(n.graph, first_bh, nv(n.graph))
    filter!(i -> i ≠ sink(n), boreholes)
    boreholes
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

"""
    Valve{T <: Number}(
        split::Dict{Int, T}
    )

Models the behaviour of the fluid at pipe divergences. `split` is a dicitonary representing the divergence.
Its keys represent the boreholes where the flow splits and its values are the proportion of flow going to each borehole. 
"""
struct Valve{T <: Number}
    split::Dict{Int, T}
end
act(v::Valve, i, m) = m * v.split[i]

"""
    equal_valve(nodes::Vector{Int}) 

Return a `Valve` that splits the flow equally to all the boreholes `nodes` in the divergence. 
"""
function equal_valve(nodes::Vector{Int}) 
    valve = Valve(Dict{Int, Float64}())
    for n in nodes
        valve.split[n] = 1/length(nodes)
    end
    return valve
end
"""
    valve(nodes::Vector{Int}, weights::Vector{T}) 

Return a `Valve` that splits the flow to each borehole `nodes[i]` in the proportion `weights[i]`.
"""
function valve(nodes::Vector{Int}, weights::Vector{T}) where {T <: Number}
    valve = Valve(Dict{Int, T}())
    for i in 1:length(nodes)
        valve.split[nodes[i]] = weights[i] 
    end
    return valve
end

"""
    absolute_valve(nodes::Vector{Int}, mass_flows::Vector{T}) 

Return a `Valve` that splits the flow such that each borehole `nodes[i]` receives an absolute amount 
of mass flow equal to `mass_flows[i]`.
"""
function absolute_valve(nodes::Vector{Int}, mass_flows::Vector{T}) where {T <: Number}
    valve = Valve(Dict{Int, T}())
    total_mass_flow = sum(mass_flows)
    for (i, n) in enumerate(nodes)
        valve.split[n] = total_mass_flow != 0 ? mass_flows[i] / total_mass_flow : 0.
    end
    return valve
end

"""
    BoreholeOperation{T <: Number}(
        valves::Dict{Int, Valve{T}}
        source_mass_flow::T
        network::BoreholeNetwork
    )

Representation of the a state of operation of the borefield. 

# Arguments
- `valves`: A dicitonary specifying the divergences of the network. Each key is the borehole 
    where the divergence occurs, and its corresponding value is a `Valve` object, describing the behaviour
    of the divergence.
- `source_mass_flow`: The total amount of mass flow injected into all branches.
- `network`: The representation of the pipe layout of the borefield, including direction of flow.

# Convenience constructors
`BoreholeOperation(;network::BoreholeNetwork, mass_flows::Vector{T})`
Builds the operation corresponding to `network` sending mass flow `mass_flows[i]` in each branch `identical`. 
This method takes care of the initial valve from the source to the start of each branch.
"""
struct BoreholeOperation{T <: Number}
    valves::Dict{Int, Valve{T}}
    source_mass_flow::T
    network::BoreholeNetwork
end

BoreholeOperation(::Nothing) = BoreholeOperation(Dict{Int, Valve{Float64}}(), 0., BoreholeNetwork())
function BoreholeOperation(;network::BoreholeNetwork, mass_flows::Vector{T}) where {T <: Number}
    source_valve = absolute_valve(outneighbors(network.graph, source(network)), mass_flows)
    BoreholeOperation(Dict(source(network) => source_valve), sum(mass_flows), network)
end

"""
    Base.reverse(network::BoreholeNetwork) 

Create the reversed flow network from the original `network`.
"""
function Base.reverse(network::BoreholeNetwork) 
    inlets = outneighbors(network.graph, source(network))
    outlets = inneighbors(network.graph, sink(network))
    reverse_graph = Base.reverse(network.graph)
    for inlet in inlets 
        rem_edge!(reverse_graph, inlet, source(network))
    end
    for outlet in outlets 
        rem_edge!(reverse_graph, sink(network), outlet)
    end
    for inlet in inlets 
        add_edge!(reverse_graph, inlet, sink(network))
    end
    for outlet in outlets 
        add_edge!(reverse_graph, source(network), outlet)
    end
    BoreholeNetwork(reverse_graph)
end


"""
    compute_mass_flows!(mass_flows, network::BoreholeNetwork, operation::BoreholeOperation{T})

Determine the mass flow each borehole receives, as determined by the network topology `network` 
and the operation `operation` and write the result in `mass_flows`.
"""
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

