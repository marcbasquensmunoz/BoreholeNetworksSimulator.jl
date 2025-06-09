using BoreholeNetworksSimulator
using BNSPlots
using Parameters


### Define input data - part that needs to be changed by the user
# Time step in seconds
Δt = 3600. # [s]
# Number of time steps (20 years with hourly resolution in this example)
Nt = 1*8760 # [-]

# Borehole buried depth (depth at which the heat extraction starts)
D = 0. # [m]

# Borehole length (length of the part of the borehole that exchanges heat)
H = 150. #[m]

# Inlet temperature to the borehole(s) to initialize the model. Can be overwritten at each time step.
Tin = 10. # [degC]

fluid = Water()

# Number of boreholes
Nb = 2

# Ground thermal diffusivity
α = 1e-6 # [m2/s]
# Ground thermal conductivity
λ = 3. # [W/(mK)]
# Undisturbed (initial) ground temperature 
T0 = 9. # [degC]

# Define the positions of each borehole (x,y)
σ = 5.
positions = [(0., 0.), (0., σ)]


### Initialize problem - no need for user intervention
# In network_1 only borehole 1 operates, and borehole 2 does not exist/operate
network_1 = BoreholeNetwork(Nb)
connect_to_source!(network_1, 1)
connect_to_sink!(network_1, 1)
connect!(network_1, 2, 2)

# In network_2 borehole 1 and borehole 2 operate in parallel
network_2 = BoreholeNetwork(Nb)
connect_to_source!(network_2, 1)
connect_to_source!(network_2, 2)
connect_to_sink!(network_2, 1)
connect_to_sink!(network_2, 2)

configurations = [network_1, network_2]

total_mass_flow = 1. #[kg/s] or [l/s]?

method = NonHistoryMethod()
medium = GroundMedium(λ=λ, α=α, T0=T0)

# Create the borehole object 
borehole = SingleUPipeBorehole(H=H, D=D)
# Create the borefield object 
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)
# Define the boundary condition
constraint = uniform_InletTempConstraint(Tin .* ones(Nt), Nb)

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

### Run the simulation
function BoreholeNetworksSimulator.operate(op::StepOperator, step, options, X)
    @unpack mass_flows, activation_step, mass_flow_containers = op
    # Tin_step = 5.
    # options.constraint.T_in[:, step] .= Tin_step
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
# fig = monitor(containers, [1, 2], options.t, Δt = :year, display=[:Tfin, :Tfout, :Tb, :q, :mf])
fig = monitor(containers, [1, 2], options.t, Δt = :year, display=[:Tfout])


# BNSPlots.get_Tfout(containers)
print(containers.X)