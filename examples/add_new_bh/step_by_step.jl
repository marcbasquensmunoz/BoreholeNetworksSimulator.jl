"""
This script simulates the operation of a geothermal borehole system in which a single borehole
operates alone for the first 3 years, after which a second borehole is connected in parallel. 
The two boreholes then continue to operate for four more years. Both boreholes share the same
geometric characteristics.

The script can be generalized to handle a different numbers of initial and added boreholes
by adjusting the relevant parameters and network definitions.

Although this example uses equal mass flow rates for all boreholes, the script supports assigning
independent mass flow rates to each one. Similarly, while the inlet temperature is set as constant
throughout the simulation, the structure allows for implementation of time-varying inlet temperature
profiles.

Such simulations are particularly useful for scenarios where an existing system is extended or
retrofitted — e.g., upgrading a residential heat pump system by adding more boreholes,
or scaling up a prototype borehole thermal storage system.
"""

using BoreholeNetworksSimulator
using BNSPlots
using Parameters
using Statistics

# --- User defined input ---
Δt = 3600.         # Time steps in seconds (1 hour)
Nt = 7 * 8760        # Total number of time steps (7 years of hourly data)

Nt_BH2 = 3*8760    # Time steps at which borehole starts operating (after 3 years)

D = 0.              # Borehole buried depth (m)
H = 150.            # Borehole active length (m)

Nb = 2              # Number of boreholes. Obs! If you change this value you have to change the positions and networks 

α = 1e-6            # Ground thermal diffusivity (m2/s)
λ = 3.              # Ground thermal conductivity (W/m/K)
T0 = 9.             # Undisturbed ground temperature (°C)

σ = 5.                              # Distance between boreholes
positions = [(0., 0.), (0., σ)]     # Coordinates of the two boreholes

Tin = 10.           # Fluid inlet temperature (°C)
m = [0.6, 0.6]      # Fluid mass flow rates through each borehole #[kg/s]

fluid = Water()

# --- Initializing the simulation  ---

# --- Borehole network configurations ---
network_1 = BoreholeNetwork(Nb)
connect_to_source!(network_1, 1)
connect_to_sink!(network_1, 1)
connect!(network_1, 2, 2)

# Configuration 2: Both boreholes are active and connected in parallel
network_2 = BoreholeNetwork(Nb)
connect_to_source!(network_2, 1)
connect_to_source!(network_2, 2)
connect_to_sink!(network_2, 1)
connect_to_sink!(network_2, 2)

# Combine both configurations for switching later during simulation
configurations = [network_1, network_2]

# --- Define simulation components ---
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

@with_kw mutable struct StepOperator{T <: Number}
    mass_flow_containers::Vector{T} = zeros(T, 2)
    mass_flows::Vector{T}
    activation_step::Int
    Tin::T = 0.
end

update(operator::StepOperator, Tin) = operator.Tin = Tin

function BoreholeNetworksSimulator.operate(op::StepOperator, step, options, X)
    @unpack mass_flows, activation_step, mass_flow_containers = op
    after_step = step >= activation_step
    active_configuration = after_step ? 2 : 1
    active_network = options.configurations[active_configuration]

    options.constraint.T_in[:, step] .= op.Tin

    if after_step 
        mass_flow_containers .= mass_flows
    else 
        mass_flow_containers[1] = mass_flows[1]
        mass_flow_containers[2] = 0.
    end
    BoreholeOperation(network=active_network, mass_flows=mass_flow_containers)
end
operator = StepOperator{Float64}(mass_flows = m, activation_step = Nt_BH2, Tin=Tin)

reset!(options)

# --- Run simulation ---
for i = 1:Nt
    if i > 1
        if i < Nt_BH2
            @views Tout = mean(containers.X[2, i-1])
        else
            @views Tout = mean(containers.X[2:2:4, i-1])
        end
        Tin = Tout - 3.
        update(operator, Tin)
    end
    simulate_steps!(n = 1, initial_step = i, operator=operator, options=options, containers=containers)
end

# --- Plot results ---
fig = monitor(containers, [1, 2], options.t, Δt = :year, display=[:Tfin, :Tfout, :Tb, :q, :mf])