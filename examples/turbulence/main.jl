using BoreholeNetworksSimulator
using BNSPlots
using Graphs
using WGLMakie

Δt = 3600.
Nt = 8760

D = 0.
H = 100.

Nb = 1

α = 1e-6
λ = 3.

laminar_mf = 0.093
turbulent_mf = 0.094

network = all_parallel_network(Nb)
positions = [(0., 0.)]
configurations = [network]

method = NonHistoryMethod()
medium = GroundMedium(λ=λ, α=α, T0=9.)
borehole = SingleUPipeBorehole(H=H, D=D, pipe_position = ((0.03, 0.0), (-0.03, 0.0)))
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)
constraint = constant_InletTempConstraint(10 .* ones(Nb), Nt)
fluid = Water()

options = SimulationOptions(
    method = method,
    constraint = constraint,
    borefield = borefield,
    fluid = fluid,
    medium = medium,
    boundary_condition = DirichletBoundaryCondition(),
    Δt = Δt,
    Nt = Nt,
    configurations = configurations
)

struct VariableMassFlowOperator <: Operator end

function BoreholeNetworksSimulator.operate(::VariableMassFlowOperator, step, options, X)
    mass_flow = floor((step / (24*30)))%2 == 0 ? laminar_mf : turbulent_mf
    network = options.configurations[1]
    source_valve = absolute_valve(Graphs.outneighbors(network.graph, source(network)), [mass_flow])
    BoreholeOperation(Dict(source(network) => source_valve), mass_flow, network)
end

operator = VariableMassFlowOperator()
containers = @time initialize(options)
@time simulate!(operator=operator, options=options, containers=containers)

fig = monitor(containers, [1], options.t, display = [:Tfout, :Tb, :q, :mf])

save("examples/turbulence/turbulence.png", fig)
