using BoreholeNetworksSimulator
using BNSPlots
using Parameters

Δt = 3 * 60.
Nt = 20*8760*20

D = 0.
H = 150.

Nb = 2

α = 1e-6
λ = 3.
T0 = 9.

σ = 5.
positions = [(0., 0.), (0., σ)]

network_1 = BoreholeNetwork(2)
connect_to_source!(network_1, 1)
connect_to_sink!(network_1, 1)
connect!(network_1, 2, 2)

network_2 = BoreholeNetwork(2)
connect_to_source!(network_2, 1)
connect_to_source!(network_2, 2)
connect_to_sink!(network_2, 1)
connect_to_sink!(network_2, 2)

configurations = [network_1, network_2]

total_mass_flow = 1.

method = NonHistoryMethod()
medium = GroundMedium(λ=λ, α=α, T0=T0)
borehole = SingleUPipeBorehole(H=H, D=D)
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)
constraint = uniform_InletTempConstraint(10. .* ones(Nt), Nb)
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

@with_kw struct StepOperator{T <: Number}
    mass_flow_containers::Vector{T} = zeros(T, 2)
    mass_flows::Vector{T}
    activation_step::Int
end

function BoreholeNetworksSimulator.operate(op::StepOperator, step, options, X)
    @unpack mass_flows, activation_step, mass_flow_containers = op
    after_step = step >= activation_step
    active_configuration = after_step ? 2 : 1
    active_network = options.configurations[active_configuration]
    if after_step 
        mass_flow_containers .= mass_flows
    else 
        mass_flow_containers[1] = mass_flows[1]
        mass_flow_containers[2] = 0.
    end
    BoreholeOperation(network=active_network, mass_flows=mass_flow_containers)
end
operator = StepOperator{Float64}(mass_flows = [0.3, 0.3], activation_step = 8760*10)

reset!(options)
@time simulate!(operator=operator, options=options, containers=containers)
fig = monitor(containers, [1, 2], options.t, Δt = :year, display=[:Tfin, :Tfout, :Tb, :q, :mf])
