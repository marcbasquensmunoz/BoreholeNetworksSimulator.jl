import sys, os
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-3]))
import BNSPythonAdapter.src.adapter
import numpy as np
import pandas as pd

borehole_positions_file = os.path.join(os.getcwd(), "Braedstrup_borehole_coordinates.txt")
Δt = 8760*3600/12
Nt = 120

network = jl.BoreholeNetwork([
    [22,30,37,29,36,35], 
    [34,40,41,42,48,43],  
    [33,32,39,45,46,47],                        
    [26,25,24,31,38,44],  
    [23,16,17,11,12,18],  
    [10,5,6,7,13,19],     
    [1,2,3,8,14,20],      
    [4,9,15,21,28,27] 
])
configurations = jl.Array([
    network,
    jl.reverse(network)
])

with open(borehole_positions_file) as csv_file:
    csv_reader = pd.read_csv(csv_file, delimiter=';', decimal=".")
    df = pd.DataFrame(csv_reader) 

borehole_positions = jl.Array[jl.Tuple[jl.Float64, jl.Float64]]([(x, y) for (x,y) in zip(df.X,df.Y)])
borehole = jl.SingleUPipeBorehole(H=50., D=4., λg = 2.5, pipe_position = ((0.03, 0.0), (-0.03, 0.0)))
borefield = jl.EqualBoreholesBorefield(borehole_prototype=borehole, positions=borehole_positions)
medium = jl.FlowInPorousMedium(λw = 0.6, λs = 2., Cw = 4.18*1e6, Cs = 1.7*1e6, θ = 0., Φ = 0.2, T0 = 10.)
constraint = jl.uniform_InletTempConstraint(jl.Array[jl.Float64]([90. if i%12 in range(1,7) else 55. for i in range(1, Nt+1)]), jl.BoreholeNetworksSimulator.n_branches(network))
method = jl.ConvolutionMethod()
fluid = jl.Fluid(cpf = 4182., name = "INCOMP::MEA-20%")

options = jl.SimulationOptions(
    method = method,
    constraint = constraint,
    borefield = borefield,
    fluid = fluid,
    medium = medium,
    Δt = Δt,
    Nt = Nt,
    configurations = configurations
)


class SeasonalOperator(jl.Operator):
    def __init__(self, mass_flows, seasonal_configuration):
        self.mass_flows = mass_flows
        self.seasonal_configuration = seasonal_configuration

    def operate(self, i, options, Tfin, Tfout, Tb, q):
        active_network = options.configurations[operator.seasonal_configuration[i]]
        jl.BoreholeOperation(active_network, operator.mass_flow)


#operator = SeasonalOperator(mass_flows=0.5 .* ones(n_branches(network)), seasonal_configuration=[i%12 in 1:6 ? 2 : 1 for i in 1:Nt])

def operator(i, Tin, Tout, Tb, q, configurations):
    mf = 0.5
    active_network = configurations[0 if i%12 in range(1, 7) else 1]
    Nbr = jl.BoreholeNetworksSimulator.n_branches(active_network)
    op = jl.BoreholeOperation(network=active_network, mass_flows=jl.Array[jl.Float64](mf * np.ones(Nbr)))
    return op

containers = jl.initialize(options)
jl.simulate_b(operator=operator, options=options, containers=containers)

containers.X