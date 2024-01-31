from juliacall import Main as jl
from juliacall import Pkg as jlPkg

jlPkg.activate("../../")
jl.seval("using BTESGroundWaterSimulator")
jl.seval("using GeometryTypes")

import numpy as np
import os
import pandas as pd

networks = jl.Array(
[
    [   
        [22,30,37,29,36,35], 
        [34,40,41,42,48,43],  
        [33,32,39,45,46,47],                        
        [26,25,24,31,38,44],  
        [23,16,17,11,12,18],  
        [10,5,6,7,13,19],     
        [1,2,3,8,14,20],      
        [4,9,15,21,28,27]                                         
    ]
    ,
    [
        [35,36,29,37,30,22],
        [43,48,42,41,40,34],
        [47,46,45,39,32,33],
        [44,38,31,24,25,26],
        [18,12,11,17,16,23],
        [19,13,7,6,5,10],
        [20,14,8,3,2,1],
        [27,28,21,15,9,4]
    ]
]
)

tstep = 8760*3600/12.
tmax  = 8760*3600*10.
Nt = int(tmax // tstep)

def operator(i, Tin, Tout, Tb, Δq, Q):
    return jl.BoreholeOperation(network=jl.Array(networks[1 if i%12 in range(6) else 0]), mass_flows=jl.Array[jl.Float64](0.5 * np.ones(8)), cpf=4182.)

borehole_positions_file = os.path.join(os.getcwd(), "../example1/data/Braedstrup_borehole_coordinates.txt")

with open(borehole_positions_file) as csv_file:
    csv_reader = pd.read_csv(csv_file, delimiter=';', decimal=",")
    df = pd.DataFrame(csv_reader) 

borehole_positions = jl.Array([jl.Point2(x, y) for (x,y) in zip(df.X,df.Y)])
borefield = jl.EqualBoreholesBorefield(borehole_prototype=jl.SingleUPipeBorehole(H=50., D=4.), positions=borehole_positions, medium=jl.GroundWaterMedium())


cache = ""

parameters = jl.compute_parameters(borefield=borefield, tstep=tstep, tmax=tmax)
constraint = jl.InletTempConstraint(jl.Array[jl.Float64]([90. if i%12 in range(6) else 55. for i in range(Nt)]))
method = jl.ConvolutionMethod(T0 = 10., parameters=parameters, borefield=borefield)
containers = jl.SimulationContainers(parameters)

jl.simulate(parameters=parameters, containers=containers, operator=operator, borefield=borefield, constraint=constraint, method=method)