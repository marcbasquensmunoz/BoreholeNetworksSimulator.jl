from juliacall import Main as jl
from juliacall import Pkg as jlPkg

jlPkg.activate("../..")
jl.seval("using BTESGroundWaterSimulator")

borefield = jl.BorefieldProperties(1. * 1e-2, 0.6, 2., 4.18*1e6, 1.7*1e6, 0., 0.2)
borehole = jl.BoreholeProperties(50., 50., 4., 2.5, 0.5, 1000., 4182.)

external_to_internal = jl.Array([ 
                        [22,30,37,29,36,35], 
                        [34,40,41,42,48,43],  
                        [33,32,39,45,46,47],                        
                        [26,25,24,31,38,44],  
                        [23,16,17,11,12,18],  
                        [10,5,6,7,13,19],     
                        [1,2,3,8,14,20],      
                        [4,9,15,21,28,27]                                         
                        ])

jl.sim1(borefield, borehole, external_to_internal, tstep = 8760*3600/12., tmax = 8760*3600*10.)
