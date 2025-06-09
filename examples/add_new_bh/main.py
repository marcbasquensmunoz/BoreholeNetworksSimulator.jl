import sys, os
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-3]))
import BNSPythonAdapter.src.adapter
from juliacall import Main as jl
import numpy as np
import pandas as pd



### Define input data - part that needs to be changed by the user
# Time step in seconds
dt = 3600. # [s]
# Number of time steps (1 year with hourly resolution in this example)
Nt = 15*8760 # [-]

# Borehole buried depth (depth at which the heat extraction starts)
D = 0. # [m]

# Borehole length (length of the part of the borehole that exchanges heat)
H = 150. #[m]

# Inlet temperature to the borehole(s) to initialize the model. Can be overwritten at each time step.
Tin = 10. # [degC]

fluid = jl.Water()

# Number of boreholes
Nb = 2

# Ground thermal diffusivity
alpha = 1e-6 # [m2/s]
# Ground thermal conductivity
λ = 3. # [W/(mK)]
# Undisturbed (initial) ground temperature 
T0 = 9. # [degC]

# Define the positions of each borehole (x,y)
σ = 5.

positions = jl.Array[jl.Tuple[jl.Float64, jl.Float64]]([(0., 0.), (0., σ)])

# total_mass_flow = 1. #[kg/s] or [l/s]?
mass_flows = jl.Array[jl.Float64](0.3 * np.ones(Nb))
activation_step = 8760*10

### Initialize problem - no need for user intervention
# In network_1 only borehole 1 operates, and borehole 2 does not exist/operate
network_1 = jl.BoreholeNetwork(Nb)
jl.connect_to_source_b(network_1, 1)
jl.connect_to_sink_b(network_1, 1)
jl.connect_b(network_1, 2, 2)

# In network_2 borehole 1 and borehole 2 operate in parallel
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
# print(constraint)

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

    def operate(self,step, options,X):
        after_step = step >= self.activation_step
        active_configuration = 1 if after_step else 0
        active_network = options.configurations[active_configuration]

        if after_step:
            self.mass_flow_containers[:] = self.mass_flows  
        else:
            self.mass_flow_containers[0] = self.mass_flows[0]  
            self.mass_flow_containers[1] = 0.0

        return jl.BoreholeOperation(network=active_network, mass_flows=operator.mass_flow_containers)


operator = StepOperator(mass_flows, activation_step,Nb)

jl.simulate_b(operator=operator, options=options, containers=containers)

# Print outlet temperature
# print(containers.X[2 * Nb, :])


### Postprocessing
import plotly.express as px

fig=px.scatter(np.array(containers.X[2 * Nb, :]))
fig.show()