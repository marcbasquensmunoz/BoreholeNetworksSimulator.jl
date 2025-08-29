using BoreholeNetworksSimulator
using BNSPlots
using Graphs
using WGLMakie

Δt = 3600*24.
Nt = 10*365

D = 0.
H = 100.

Nb = 2

α = 1e-6
λ = 3.
T0 = 9.

σ = 1.

network = all_series_network(Nb)
positions = [(0., 0.), (0., σ)]
configurations = [network, reverse(network)]

Q = H*Nb
total_mass_flow = 1.

method = OriginalNonHistoryMethod()
medium = GroundMedium(λ=λ, α=α, T0=T0)
borehole = SingleUPipeBorehole(H=H, D=D)
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)
constraint = TotalHeatLoadConstraint(Q * ones(Nt))
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
containers = @time initialize(options)

struct ReversalStrategy
    mass_flows
    Q_tot
end

function BoreholeNetworksSimulator.operate(op::ReversalStrategy, step, options, X)
    @unpack Q_tot, mass_flows = op
    active_configuration = floor((step + 3*30) / (6*30)) % 2 == 0 ? 1 : 2
    injection = floor((step) / (6*30)) % 2 == 0 ? -1 : 1
    options.constraint.Q_tot .= injection * Q_tot
    active_network = options.configurations[active_configuration]
    BoreholeOperation(network=active_network, mass_flows=mass_flows)
end

operator = ReversalStrategy(total_mass_flow/2 .* ones(Nb), Q)

reset!(options)
@time simulate!(operator=operator, options=options, containers=containers)
fig = monitor(containers, [1, 2], options.t, Δt = :year, display=[:Tb, :q])
save("$(@__DIR__)/plots/reversal.png", fig)
