using BoreholeNetworksSimulator
using BNSPlots
using Graphs
using WGLMakie

Δt = 3600.
Nt = 8760

D = 0.
H = 100.

Nb = 2

α = 1e-6
λ = 3.
T0 = 9.

σ = 5.

network = all_parallel_network(Nb)
positions = [(0., 0.), (0., σ)]
configurations = [network]

Q = H
total_mass_flow = 1.

method = NonHistoryMethod()
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
