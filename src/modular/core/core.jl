
@with_kw struct BoreholeOperation
    network         
    mass_flows:: Vector          # Mass flow for each branch in network
end
BoreholeOperation(::Nothing) = BoreholeOperation([[]], [0.])

@with_kw struct Fluid
    cpf
    name
end

@with_kw struct BoreholeNetwork
    branches::Vector
end
Base.reverse(network::BoreholeNetwork) = BoreholeNetwork(branches=map(branch -> Base.reverse(branch), network.branches))
n_branches(network::BoreholeNetwork) = length(network.branches)
first_boreholes(network::BoreholeNetwork) = map(first, network.branches)

@with_kw struct SimulationOptions
    method::TimeSuperpositionMethod
    constraint::Constraint
    borefield::Borefield
    fluid::Fluid
    boundary_condition::BoundaryCondition = DirichletBoundaryCondition()
    Δt
    Nt
    Nb = n_boreholes(borefield)
    Ns = n_segments(borefield)
    Ts = 1
    Tmax = Δt * Nt
    t = Δt:Δt:Tmax
end

@with_kw struct SimulationContainers
    M
    X
    b
end

function initialize(options::SimulationOptions) 
    precompute_auxiliaries!(options.method, options=options)
    SimulationContainers(options)
end
function SimulationContainers(options::SimulationOptions) 
    @unpack Nb, Nt = options
    SimulationContainers(M = spzeros(4Nb, 4Nb), b = zeros(4Nb), X = zeros(4Nb, Nt))
end

function branch_of_borehole(operation::BoreholeOperation, borehole)
    for (i, branch) in enumerate(operation.network.branches)
        if borehole in branch
            return i
        end
    end
    return 0
end


function solve_step!(X, A, b)
    #=
    prob = LinearProblem(A, b)
    linsolve = init(prob)
    x = solve(linsolve).u
    =#
    X .= A\b
end

function topology_coeffs!(M, operation::BoreholeOperation)
    for (i, branch) in enumerate(operation.network.branches)
        for (out, in) in zip(branch[1:end-1], branch[2:end])
            M[i, 2*in-1] = 1.
            M[i, 2*out] = -1.
        end
    end
end
