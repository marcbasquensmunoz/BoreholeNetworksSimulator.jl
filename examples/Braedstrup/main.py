import sys, os
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-3]))
import BNSPythonAdapter.src.adapter
from juliacall import Main as jl
import numpy as np
import pandas as pd

borehole_positions_file = os.path.join(os.getcwd(), "data", "Braedstrup_borehole_coordinates.txt")
Δt = 8760*3600/12
Nt = 120

network = jl.BoreholeNetwork(48)
jl.connect_to_source_b(network, jl.Vector[jl.Int]([35, 43, 47, 44, 18, 19, 20, 27]))
jl.connect_in_series_b(network, jl.Vector[jl.Int]([35, 36, 29, 37, 30, 22]))
jl.connect_in_series_b(network, jl.Vector[jl.Int]([43, 48, 42, 41, 40, 34]))
jl.connect_in_series_b(network, jl.Vector[jl.Int]([47, 46, 45, 39, 32, 33]))
jl.connect_in_series_b(network, jl.Vector[jl.Int]([44, 38, 31, 24, 25, 26]))
jl.connect_in_series_b(network, jl.Vector[jl.Int]([18, 12, 11, 17, 16, 23]))
jl.connect_in_series_b(network, jl.Vector[jl.Int]([19, 13, 7, 6, 5, 10]))
jl.connect_in_series_b(network, jl.Vector[jl.Int]([20, 14, 8, 3, 2, 1]))
jl.connect_in_series_b(network, jl.Vector[jl.Int]([27, 28, 21, 15, 9, 4]))
jl.connect_to_sink_b(network, jl.Vector[jl.Int]([22, 34, 33, 26, 23, 10, 1, 4]))

configurations = jl.Vector([
    network,              # Heat extraction
    jl.reverse(network)   # Heat injection
])

with open(borehole_positions_file) as csv_file:
    csv_reader = pd.read_csv(csv_file, delimiter=';', decimal=",")
    df = pd.DataFrame(csv_reader) 


Tf_injection = 90.
Tf_extraction = 55.

borehole_positions = jl.Array[jl.Tuple[jl.Float64, jl.Float64]]([(x, y) for (x,y) in zip(df.X,df.Y)])
borehole = jl.SingleUPipeBorehole(H=50., D=4., λg = 1.5, pipe_position = ((0.03, 0.0), (-0.03, 0.0)))
borefield = jl.EqualBoreholesBorefield(borehole_prototype=borehole, positions=borehole_positions)
medium = jl.FlowInPorousMedium(λw = 0.6, λs = 2., Cw = 4.18*1e6, Cs = 1.7*1e6, θ = 0., Φ = 0.2, T0 = 10.)
constraint = jl.uniform_InletTempConstraint(jl.Array[jl.Float64]([Tf_injection if i%12 in range(1,7) else Tf_extraction for i in range(1, Nt+1)]), jl.BoreholeNetworksSimulator.n_branches(network))
method = jl.ConvolutionMethod()
fluid = jl.Water()

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


class SeasonalOperator():
    def __init__(self, mass_flows, seasonal_configuration):
        self.mass_flows = mass_flows
        self.seasonal_configuration = seasonal_configuration

    def operate(self, i, options, X):
        active_network = options.configurations[operator.seasonal_configuration[i-1]]
        return jl.BoreholeOperation(network=active_network, mass_flows=operator.mass_flows)


seasonal_configuration = [1 if i%12 in range(1, 7) else 0 for i in range(1,Nt+1)]
mass_flows = jl.Array[jl.Float64](0.5 * np.ones(jl.n_branches(network)))
operator = SeasonalOperator(mass_flows, seasonal_configuration)

containers = jl.initialize(options)
jl.simulate_b(operator=operator, options=options, containers=containers)

containers.X