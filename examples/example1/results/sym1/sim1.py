from juliacall import Main as jl
from juliacall import Pkg as jlPkg

jlPkg.activate("../BTESGroundWaterSimulator")
jl.seval("using BTESGroundWaterSimulator")

jl.seval("using BoreholeResponseFunctions")
jl.seval("using GeometryTypes")
jl.seval("using Parameters")
jl.seval("using CSV")
jl.seval("using DataFrames")
jl.seval("using LinearAlgebra")

import numpy as np
import os

# 1. PROPERTIES OF BOREHOLE FIELD
ux_in_meterperday  = 1. * 1e-2      # groundwater speed along the flow coordinate
ux =  ux_in_meterperday/(3600*24)
λw = 0.6                            # water thermal conductivity
λs = 2.                             # ground thermal conductivity
Cw = 4.18*1e6                       # water thermal capacity
Cs = 1.7*1e6                        # ground thermal capacity
θ = 0.                              # angle of Darcy velocity
Φ = 0.2                             # porosity
λ = λs *(1-Φ) + λw*Φ                # porous medium conductivity
C = Cs *(1-Φ) + Cw*Φ                # porous medium capacity
α  = λ/C                            # porus medium thermal diffusivity
vt = ux * Cw/C                      # porous medium darcy velocity

# 2. GEOMETRICAL CONFIGURATION OF THE BOREHOLES ALONG THE VERTICAL DIRECTION

# Geometry on the vertical direction 

# H length of the borehole 
# h length of segment 
# D groundwater level
H,h,D = 50., 50., 4.


# 3. U-PIPE MODEL 
params = jl.BoreholePara(λg=2.5,λs = λ)                          # Define borehole u-pipe cross-section geometry 
upipe_crossection = jl.broadcast(jl.GeometryTypes.Point2, ([[0.03,0.0],[-0.03,.0]]))   # (x,y) pipe positions on the borehole cross-section
rb  = params.rb
rpo = params.rpo

mf  =  0.5          # mass flow kg/s
ρf  =  1000.        # density 
cpf =  4182.        # specific heat capacity
Cf  =  ρf*cpf       # capacity  
Vf  =  mf/ρf        # volume flow

R = jl.resistance_network(params, upipe_crossection)   # resistance network in the cross-section
A = jl.coefficient_matrix(R, Cf, Vf)


# 4. BOREHOLE FIELD CONFIGURATION
#import configuration
cdir = os.getcwd()
cdir = os.path.join(cdir,"examples/example1")
tdir = os.path.join(cdir,"data/Braedstrup_borehole_coordinates.txt")
df = jl.DataFrame(jl.CSV.File(tdir)) 
#geometry of the field
borehole_positions =  [(x, y) for (x,y) in zip(df.X,df.Y) ]

# borehole connections map for extraction mode
external_to_internal = jl.Array([ [22,30,37,29,36,35], 
                        [34,40,41,42,48,43],  
                        [33,32,39,45,46,47],                        
                        [26,25,24,31,38,44],  
                        [23,16,17,11,12,18],  
                        [10,5,6,7,13,19],     
                        [1,2,3,8,14,20],      
                        [4,9,15,21,28,27]                                         
                        ]
                        )

# borehole connections map for injection mode
internal_to_external = jl.broadcast(jl.reverse, external_to_internal)

# 5. DISCRETIZATION
z_ref = jl.collect(jl.range(D, step = h, stop = D+H-h) )         # line source reference point
z_eval = jl.collect(jl.range(D+h/2, step = h, stop = D+H-h/2))   # evaluation points (evaluate at the mid point of the segment)

Nb = len(borehole_positions)     # number of boreholes
Nsb = len(z_eval)                # number of segments per borehole !
Ns  = Nb*Nsb                        # total number of segmensts

bh_map     =  jl.reshape(jl.ones(jl.Int64, Nsb)* jl.transpose(jl.collect(jl.range(1,Nb))),(Ns,1))
coord_source = [(x[0],x[1],p) for x in borehole_positions for p in z_ref] # position of sources 
coord_eval = [(x[0],x[1],p) for x in borehole_positions for p in z_eval]  # position of evaluation points

# p = jl.broadcast(jl.GeometryTypes.Point3{jl.Float64} , coord_source)
# tp = GeometryTypes.Point3{Float64}.(coord_eval)  

# # # rotation of points in new coordinate system where Darcy velocity is parallel to x axis
# p_rot  = rotation_z(p,-θ) 
# tp_rot = rotation_z(tp,-θ) 

# d = evaluate_relevant_distances(GroundWaterFlow(), p_rot, tp_rot) 
# d = [d[1] == 0. && d[2] == 0. ?  (0.,params.rb, d[3],d[4]) : d for d in d]

# tstep, tmax = 8760*3600/12., 8760*3600*10.
# t = tstep:tstep:tmax
# Nt = length(t) # number of time steps
