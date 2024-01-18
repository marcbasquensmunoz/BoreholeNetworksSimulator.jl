from juliacall import Main as jl
from juliacall import Pkg as jlPkg

jlPkg.activate("../..")
jl.seval("using BTESGroundWaterSimulator")

borefield = jl.BorefieldProperties(1. * 1e-2, 0.6, 2., 4.18*1e6, 1.7*1e6, 0., 0.2)
borehole = jl.BoreholeProperties(50., 50., 4., 2.5, 0.5, 1000., 4182.)
configurations = jl.Array([
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
])

def operator(i, X, q):
    return 2 if (i-1)%12 in range(0,6) else 1


jl.sim1(operator=operator, borefield=borefield, borehole=borehole, configurations=configurations, tstep = 8760*3600/12., tmax = 8760*3600*10.)
