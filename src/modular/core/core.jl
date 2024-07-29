
"""
    BoreholeNetwork(branches::Vector)

Representation of the hydraulic connections of the boreholes in the network.
Each element in `branches` should be a vector representing a branch of boreholes connected in series, specified by their identifiers.
The first borehole of each branch is assumed to be connected in parallel. 
"""
@with_kw struct BoreholeNetwork
    branches::Vector
end
Base.reverse(network::BoreholeNetwork) = BoreholeNetwork(branches=map(branch -> Base.reverse(branch), network.branches))
n_branches(network::BoreholeNetwork) = length(network.branches)
first_boreholes(network::BoreholeNetwork) = map(first, network.branches)

"""
    BoreholeOperation{T <: Number}(network::BoreholeNetwork, mass_flows:: Vector{T})

Represents a operation state of the network, with `network` representing the hydraulic configuration and `mass_flows` a `Vector` containing the mass flow rate of each branch.
"""
@with_kw struct BoreholeOperation{T <: Number}
    network::BoreholeNetwork         
    mass_flows::Vector{T}
end
BoreholeOperation(::Nothing) = BoreholeOperation(BoreholeNetwork([]), [0.])

"""
    Fluid(cpf, name)

Represents the fluid flowing through the hydraulic circuit.
`cpf` is the specific heat of the fluid and `name` is the code used in CoolProp
"""
@with_kw struct Fluid
    cpf
    name
end

"""
    struct SimulationOptions{
                    TSM <: TimeSuperpositionMethod,
                    C <: Constraint,
                    B <: Borefield, 
                    M <: Medium, 
                    BC <: BoundaryCondition
                }(
        method::TSM
        constraint::C
        borefield::B
        medium::M
        Δt, 
        Nt::Int,
        boundary_condition::BoundaryCondition = DirichletBoundaryCondition(),
        fluid::Fluid = Fluid(cpf=4182., name="INCOMP::MEA-20%")
    )

Specifies all the options for the simulation.
"""
@with_kw struct SimulationOptions{
                    TSM <: TimeSuperpositionMethod,
                    C <: Constraint,
                    B <: Borefield, 
                    M <: Medium, 
                    BC <: BoundaryCondition
                }
    method::TSM
    constraint::C
    borefield::B
    medium::M
    fluid::Fluid = Fluid(cpf=4182., name="INCOMP::MEA-20%")
    boundary_condition::BC = DirichletBoundaryCondition()
    Δt
    Nt::Int
    Nb::Int = n_boreholes(borefield)
    Ns::Int = n_segments(borefield)
    Ts::Int = 1
    Tmax = Δt * Nt
    t = Δt:Δt:Tmax
end

"""
    SimulationContainers(M, X, b)

Contains the matrix `M`, the independent vector `b` defining the problem. Both `M` and `b` change through time.
Each column of `X` contains the solution of `M X = b` for each time step of the simulation. 
"""
@with_kw struct SimulationContainers
    M
    X
    b
end

"""
    initialize(options::SimulationOptions) 

Precompute the objects of each `TimeSuperpositionMethod` that can be computed ahead of time and return the `SimulationContainers` of the required size.
"""
function initialize(options::SimulationOptions) 
    precompute_auxiliaries!(options.method, options)
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
