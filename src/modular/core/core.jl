
"""
    struct SimulationOptions{
                    N <: Number,
                    Tol <: Number,
                    TSM <: TimeSuperpositionMethod,
                    C <: Constraint,
                    B <: Borefield, 
                    M <: Medium, 
                    BC <: BoundaryCondition,
                    A <: Approximation,
                    F <: Fluid
                }(
        method::TSM
        constraint::C
        borefield::B
        medium::M
        boundary_condition::BoundaryCondition = DirichletBoundaryCondition()
        approximation::A = MeanApproximation()
        fluid::Fluid{N} = Fluid(cpf=4182., name="INCOMP::MEA-20%")
        configurations::Vector{BoreholeNetwork}
        Δt
        Nt::Int
        atol::Tol = 0.
        rtol::Tol = sqrt(eps())
    )

Specifies all the options for the simulation.

- `method`: time superposition method used to compute the response. Available options: `ConvolutionMethod`, `NonHistoryMethod`.
- `constraint`: constraint that the system must satisfy. Can be variable with time. Available options: `HeatLoadConstraint`, `InletTempConstraint`, `TotalHeatLoadConstraint`.
- `borefield`: describes the geometrical properties and the boreholes of the borefield on which the simulation will be performed. Available options: `EqualBoreholesBorefield`.
- `medium`: properties of the ground where the `borefield` is places. Available options: `GroundMedium`, `FlowInPorousMedium`.
- `boundary_condition`: boundary condition of the domain where the simulation is performed. Available options: `NoBoundary`, `DirichletBoundaryCondition`, `NeumannBoundaryCondition`.
- `approximation`: determines how the approximate value for each segment is computed. Available options: `MeanApproximation`, `MidPointApproximation`.
- `fluid`: properties of the fluid flowing through the hydraulic system.
- `configurations`: possible hydraulic topologies possible in the system, including reverse flow.
- `Δt`: time step used in the simulation.
- `Nt`: total amount of time steps of the simulation.
- `atol`: absolute tolerance used for the adaptive integration methods.
- `rtol`: relative tolerance used for the adaptive integration methods.
"""
@with_kw struct SimulationOptions{
                    N <: Number,
                    Tol <: Number,
                    TSM <: TimeSuperpositionMethod,
                    C <: Constraint,
                    B <: Borefield, 
                    M <: Medium, 
                    BC <: BoundaryCondition, 
                    A <: Approximation,
                    F <: Fluid
                }
    method::TSM
    constraint::C
    borefield::B
    medium::M
    fluid::F
    boundary_condition::BC = DirichletBoundaryCondition()
    approximation::A = MeanApproximation()
    Δt::N
    Nt::Int
    Nb::Int = n_boreholes(borefield)
    Ns::Int = n_segments(borefield)
    Ts::Int = 1
    Tmax::N = Δt * Nt
    t::Vector{N} = collect(Δt:Δt:Tmax)
    configurations::Vector{BoreholeNetwork}
    atol::Tol = 0.
    rtol::Tol = sqrt(eps())
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
    mf::Matrix{T}
end

Base.copy(c::SimulationContainers) = SimulationContainers(c.M, copy(c.X), copy(c.b), copy(c.mf))

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
    SimulationContainers(M = Nb > 100 ? spzeros(4Nb, 4Nb) : zeros(4Nb, 4Nb), b = zeros(4Nb), X = zeros(4Nb, Nt), mf = zeros(Nb, Nt))
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

function reset!(options::SimulationOptions)
    reset!(options.method)
end

reset!(::TimeSuperpositionMethod) = nothing
