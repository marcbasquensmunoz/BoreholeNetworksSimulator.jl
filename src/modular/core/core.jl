
"""
    BoreholeNetwork(branches::Vector{Vector{Int}})

Representation of the hydraulic connections of the boreholes in the network.
Each element in `branches` should be a vector representing a branch of boreholes connected in series, specified by their identifiers.
The first borehole of each branch is assumed to be connected in parallel. 
"""
@with_kw struct BoreholeNetwork
    branches::Vector{Vector{Int}}
end
Base.reverse(network::BoreholeNetwork) = BoreholeNetwork(branches=map(branch -> Base.reverse(branch), network.branches))
n_branches(network::BoreholeNetwork) = length(network.branches)
n_boreholes(network::BoreholeNetwork) = sum([length(branch) for branch in network.branches])
first_boreholes(network::BoreholeNetwork) = map(first, network.branches)

"""
    BoreholeOperation{Arr <: AbstractArray}(
        network::BoreholeNetwork         
        mass_flows::Arr
    )

Represents a operation state of the network, with `network` representing the hydraulic configuration and `mass_flows` a `Vector` containing the mass flow rate of each branch.
"""
@with_kw struct BoreholeOperation{Arr <: AbstractArray}
    network::BoreholeNetwork         
    mass_flows::Arr
end
BoreholeOperation(::Nothing) = BoreholeOperation(BoreholeNetwork([]), @view ones(1)[1:1])

"""
    Fluid{T <: Number}(
        cpf::T
        name::String
    )

Represents the fluid flowing through the hydraulic circuit.
`cpf` is the specific heat of the fluid and `name` is the code used in CoolProp
"""
@with_kw struct Fluid{T <: Number}
    cpf::T
    name::String
end

"""
    struct SimulationOptions{
                    TSM <: TimeSuperpositionMethod,
                    C <: Constraint,
                    B <: Borefield, 
                    M <: Medium, 
                    BC <: BoundaryCondition,
                    N <: Number
                }(
        method::TSM
        constraint::C
        borefield::B
        medium::M
        boundary_condition::BoundaryCondition = DirichletBoundaryCondition()
        fluid::Fluid{N} = Fluid(cpf=4182., name="INCOMP::MEA-20%")
        configurations::Vector{BoreholeNetwork}
        Δt
        Nt::Int
    )

Specifies all the options for the simulation.

- `method`: time superposition method used to compute the response. Available options: `ConvolutionMethod`, `NonHistoryMethod`.
- `constraint`: constraint that the system must satisfy. Can be variable with time. Available options: `HeatLoadConstraint`, `InletTempConstraint`.
- `borefield`: describes the geometrical properties and the boreholes of the borefield on which the simulation will be performed. Available options: `EqualBoreholesBorefield`.
- `medium`: properties of the ground where the `borefield` is places. Available options: `GroundMedium`, `FlowInPorousMedium`.
- `boundary_condition`: boundary condition of the domain where the simulation is performed. Available options: `NoBoundary`, `DirichletBoundaryCondition`.
- `fluid`: properties of the fluid flowing through the hydraulic system.
- `configurations`: possible hydraulic topologies possible in the system, including reverse flow.
- `Δt`: time step used in the simulation.
- `Nt`: total amount of time steps of the simulation.
"""
@with_kw struct SimulationOptions{
                    TSM <: TimeSuperpositionMethod,
                    C <: Constraint,
                    B <: Borefield, 
                    M <: Medium, 
                    BC <: BoundaryCondition, 
                    N <: Number
                    #A <: Approximation
                }
    method::TSM
    constraint::C
    borefield::B
    medium::M
    fluid::Fluid{N} = Fluid(cpf=4182., name="INCOMP::MEA-20%")
    boundary_condition::BC = DirichletBoundaryCondition()
    #approximation::A = MeanApproximation()
    Δt
    Nt::Int
    Nb::Int = n_boreholes(borefield)
    Ns::Int = n_segments(borefield)
    Ts::Int = 1
    Tmax = Δt * Nt
    t = Δt:Δt:Tmax
    configurations::Vector{BoreholeNetwork}
end

"""
    SimulationContainers{T <: Number, Mat <: AbstractMatrix}(
        M::Mat
        X::Matrix{T}
        b::Vector{T}
    )

Contains the matrix `M`, the independent vector `b` defining the problem. Both `M` and `b` change through time.
Each column of `X` contains the solution of `M X = b` for each time step of the simulation. 
"""
@with_kw struct SimulationContainers{T <: Number, Mat <: AbstractMatrix}
    M::Mat
    X::Matrix{T}
    b::Vector{T}
end

"""
    initialize(options::SimulationOptions) 

Precompute the objects of each `TimeSuperpositionMethod` that can be computed ahead of time and return the `SimulationContainers` of the required size.
"""
function initialize(options::SimulationOptions) 
    # TODO: detect overlapping boreholes
    compatibility = check_compatibility(options.medium, options.constraint, options.method)
    if compatibility isa NotCompatible
        println(compatibility.message)
        return
    end
    precompute_auxiliaries!(options.method, options)
    SimulationContainers(options)
end
function SimulationContainers(options::SimulationOptions) 
    @unpack Nb, Nt = options
    if Nb > 50
        SimulationContainers(M = spzeros(4Nb, 4Nb), b = zeros(4Nb), X = zeros(4Nb, Nt))
    else 
        SimulationContainers(M = zeros(4Nb, 4Nb), b = zeros(4Nb), X = zeros(4Nb, Nt))
    end
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
    X .= solve!(linsolve).u=#
    X .= A\b
end

unwrap(x::Any) = x