using Parameters
using SparseArrays

@with_kw struct SimulationParameters
    Nb
    Ns
    tstep
    tmax
    t = tstep:tstep:tmax
    Nt = length(t)
    Ts = 1
end

@with_kw struct SimulationContainers
    M
    X
    b
    current_Q
end
function SimulationContainers(parameters::SimulationParameters) 
    @unpack Nb, Ns, Nt = parameters
    SimulationContainers(M = spzeros(3Nb + Ns, 3Nb + Ns), b = zeros(3Nb + Ns), X = zeros(Nt, 3Nb + Ns), current_Q = zeros(Ns))
end
          
@with_kw struct BoreholeOperation{T <: Real}
    network::Vector{Vector{Int}}            
    mass_flows:: Vector{T}          # Mass flow for each branch in network
    cpf::T = 4182.                  # specific heat capacity
end
BoreholeOperation(::Nothing) = BoreholeOperation([[]], [0.], 0.)

function branch_of_borehole(operation::BoreholeOperation, borehole)
    for (i, branch) in enumerate(operation.network)
        if borehole in branch
            return i
        end
    end
    return 0
end

function heat_balance_coeffs!(M, borefield::Borefield, operation::BoreholeOperation)
    Nb = borehole_amount(borefield)
    Ns = segment_amount(borefield)

    for i in 1:Nb
        heat = operation.cpf .* operation.mass_flows[branch_of_borehole(operation, i)] 
        M[i, i*2-1] = heat
        M[i, i*2]  = -heat
    end

    for i in 1:Nb
        for j in 1:Ns
            if where_is_segment(borefield, j) == i
                M[i, 3Nb+j] = -get_h(borefield, i)
            end
        end
    end   
end

function heat_balance_b!(b, borefield, current_Q)
    Nb = borehole_amount(borefield)
    for i in 1:Nb
        b[i] = current_Q[i] * get_H(borefield, i)
    end  
end

function solve_step!(X, A, b, step, Nb, current_Q)
    x = A\b
    X[step,:] = x
    current_Q .+= @views x[3Nb+1:end]
end

function compute_parameters(;borefield::Borefield, tstep, tmax)
    SimulationParameters(Nb=borehole_amount(borefield), Ns=segment_amount(borefield), tstep=tstep, tmax=tmax)
end

function load_cache!(;containers::SimulationContainers, parameters::SimulationParameters, cache)
    if cache != ""
        @unpack X, b, current_Q = containers
        data = load(cache)
        parameters.Ts = size(data["X"])[1]
        X[1:parameters.Ts, :] = data["X"]
        b = data["b"]
        current_Q = data["current_Q"]
    end
end

function save_cache(;containers::SimulationContainers, parameters::SimulationParameters, path, title)
    results_directory = "$path/results"
    simulation_results_directory = "$results_directory/$title"
    !isdir(results_directory) && mkdir(results_directory)
    !isdir(simulation_results_directory) && mkdir(simulation_results_directory)
    save("$(simulation_results_directory)/cache_$(parameters.tmax).jld2" , 
        Dict( 
            "X" => containers.X,
            "b" => containers.b,
            "current_Q" => containers.current_Q
        )
    )
end