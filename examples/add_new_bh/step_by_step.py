"""
This script simulates the operation of a geothermal borehole system in which a single borehole
operates alone for the first year, after which a second borehole is connected in parallel. 
The two boreholes then continue operate for two more years. Both boreholes share the same
geometric characteristics.

The script  can be generalized to handle different numbers of initial and added boreholes
by adjusting the relevant parameters and network definitions.

Although this example uses equal mass flow rates for all boreholes, the script supports assigning
independent mass flow rates to each one. The inlet temperature is time dependent, and in particular 
it is set as equal to the outlet temperature - 3 K.

Such simulations are particularly useful for scenarios where an existing system is extended or
retrofitted — e.g., upgrading a residential heat pump system by adding more boreholes,
or scaling up a prototype borehole thermal storage system.
"""

import sys, os
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-3]))
import BNSPythonAdapter.src.adapter
from juliacall import Main as jl
import numpy as np
import pandas as pd

# --- User defined input ---
dt = 3600.          # Time steps in seconds (1 hour)
Nt = 3 * 8760         # Total number of time steps (20 years of hourly data)

Nt_BH2 = 1 * 8760     # Time steps at which borehole starts operating (after 10 years)

D = 0.              # Borehole buried depth (m)
H = 150.            # Borehole active length (m)

Nb = 2              # Number of boreholes. Obs! If you change this value you have to change the positions and networks 

alpha = 1e-6            # Ground thermal diffusivity (m2/s)
λ = 3.                  # Ground thermal conductivity (W/m/K)
T0 = 9.                 # Undisturbed ground temperature (°C)

σ = 5.                              # Distance between boreholes

positions = jl.Array[jl.Tuple[jl.Float64, jl.Float64]]([(0., 0.), (0., σ)])     # Coordinates of the two boreholes

Tin = T0 - 3.                                             # Initial fluid inlet temperature (°C)
mass_flows = jl.Array[jl.Float64](0.3 * np.ones(Nb))      # Fluid mass flow rates through each borehole #[kg/s]

fluid = jl.Water()

# --- Define two different borehole network configurations ---
# Configuration 1: Borehole 1 is active, Borehole 2 is inactive
network_1 = jl.BoreholeNetwork(Nb)
jl.connect_to_source_b(network_1, 1)
jl.connect_to_sink_b(network_1, 1)
jl.connect_b(network_1, 2, 2)

# Configuration 2: Both boreholes are active and connected in parallel
network_2 = jl.BoreholeNetwork(Nb)
jl.connect_to_source_b(network_2, 1)
jl.connect_to_source_b(network_2, 2)
jl.connect_to_sink_b(network_2, 1)
jl.connect_to_sink_b(network_2, 2)

configurations = jl.Vector([network_1, network_2])

method = jl.NonHistoryMethod()
medium = jl.GroundMedium(λ=λ, α=alpha, T0=T0)

# Create the borehole object 
borehole = jl.SingleUPipeBorehole(H=H, D=D)

# Create the borefield object 
borefield = jl.EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)
# Define the boundary condition
constraint = jl.uniform_InletTempConstraint(jl.Array[jl.Float64]([Tin for i in range(1, Nt+1)]), Nb)

options = jl.SimulationOptions(
    method = method,
    constraint = constraint,
    borefield = borefield,
    fluid = fluid,
    medium = medium,
    boundary_condition = jl.DirichletBoundaryCondition(),
    Δt = dt,
    Nt = Nt,
    configurations = configurations
)

containers = jl.initialize(options)

class StepOperator():
    def __init__(self, mass_flows, activation_step,Nb):
        self.mass_flow_containers = jl.Array[jl.Float64](np.zeros(Nb))
        self.mass_flows = mass_flows
        self.activation_step = activation_step
        self.Tin = 0

    def update(self, Tin):
        self.Tin = Tin

    def operate(self,step, options, X):
        after_step = step >= self.activation_step
        active_configuration = 1 if after_step else 0
        active_network = options.configurations[active_configuration]

        options.constraint.T_in[:, step] = self.Tin

        if after_step:
            self.mass_flow_containers[:] = self.mass_flows  
        else:
            self.mass_flow_containers[0] = self.mass_flows[0]  
            self.mass_flow_containers[1] = 0.0

        return jl.BoreholeOperation(network=active_network, mass_flows=operator.mass_flow_containers)


operator = StepOperator(mass_flows, Nt_BH2, Nb)

# --- Run simulation ---
for i in range(0,Nt):

    if i > 0:
        if i < Nt_BH2:
            Tout = np.mean(containers.X[1, i - 1])
        else: 
            Tout = np.mean(containers.X[1:4:2, i - 1])

        Tin = Tout - 3

    operator.update(Tin)
    jl.simulate_steps_b(n = 1, operator=operator, options=options, containers=containers)


# --- Plot results ---
import plotly.express as px

# Plot temperature
fig=px.scatter(np.array(containers.X[2 * Nb, :]))
fig.update_layout(yaxis=dict(range=[0, 10]))
fig.show()

# Plot inlet temperature (BH 1)
fig0=px.scatter(np.array(containers.X[0, :]))
fig0.update_layout(yaxis=dict(range=[-10, 10]))
fig0.show()

# Plot outlet temperature (BH 1)
fig1=px.scatter(np.array(containers.X[1, :]))
fig1.update_layout(yaxis=dict(range=[-10, 10]))
fig1.show()

# Plot inlet temperature (BH 2)
# When BH 2 is not operating, the borehole wall temperature is plotted instead
fig2=px.scatter(np.array(containers.X[2, :]))
fig2.update_layout(yaxis=dict(range=[-10, 10]))
fig2.show()

# Plot outlet temperature (BH 2)
# When BH 2 is not operating, the borehole wall temperature is plotted instead
fig3=px.scatter(np.array(containers.X[3, :]))
fig3.update_layout(yaxis=dict(range=[-10,10]))
fig3.show()