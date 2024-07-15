
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
end
function SimulationContainers(parameters::SimulationParameters) 
    @unpack Nb, Ns, Nt = parameters
    SimulationContainers(M = spzeros(3Nb + Ns, 3Nb + Ns), b = zeros(3Nb + Ns), X = zeros(3Nb + Ns, Nt))
end
          
@with_kw struct BoreholeOperation{T}
    network         
    mass_flows:: Vector{T}          # Mass flow for each branch in network
    cpf::T = 4182.                  # specific heat capacity
end
BoreholeOperation(::Nothing) = BoreholeOperation([[]], [0.], 0.)

@with_kw struct BoreholeNetwork
    branches::Vector
end
Base.reverse(network::BoreholeNetwork) = BoreholeNetwork(branches=map(branch -> Base.reverse(branch), network.branches))
n_branches(network::BoreholeNetwork) = length(network.branches)
first_boreholes(network::BoreholeNetwork) = map(first, network.branches)

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
    X = A\b
end

function compute_parameters(;borefield::Borefield, tstep, tmax)
    SimulationParameters(Nb=borehole_amount(borefield), Ns=segment_amount(borefield), tstep=tstep, tmax=tmax)
end

function topology_coeffs!(M, operation::BoreholeOperation)
    for (i, branch) in enumerate(operation.network.branches)
        for (out, in) in zip(branch[1:end-1], branch[2:end])
            M[i, 2*in-1] = 1.
            M[i, 2*out] = -1.
        end
    end
end
