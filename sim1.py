from juliacall import Main as jl
from juliacall import Pkg as jlPkg

jlPkg.activate(".")
jl.seval("using BTESGroundWaterSimulator")

jl.seval("using BoreholeResponseFunctions")
jl.seval("using GeometryTypes")
jl.seval("using Parameters")
jl.seval("using CSV")
jl.seval("using DataFrames")
jl.seval("using LinearAlgebra")

import numpy as np
import os
import pandas as pd

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

# determine coefficients k_in, k_out and k_b for matrix build (Cimmino (2016))
k_in, k_out, k_b = jl.uniformTb_koeff(A,H)    # COEFFICIENTS OF THE BOREHOLE MODEL


# 4. BOREHOLE FIELD CONFIGURATION
#import configuration
cdir = os.getcwd()
cdir = os.path.join(cdir,"examples/example1")
tdir = os.path.join(cdir,"data/Braedstrup_borehole_coordinates.txt")

with open(tdir) as csv_file:
    csv_reader = pd.read_csv(csv_file, delimiter=';', decimal=",")
    df = jl.DataFrame(csv_reader) 

#geometry of the field
borehole_positions =  [(x, y) for (x,y) in zip(df.X,df.Y)]

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

p = jl.broadcast(jl.GeometryTypes.Point3[jl.Float64], coord_source)
tp = jl.broadcast(jl.GeometryTypes.Point3[jl.Float64], coord_eval)  

# # # rotation of points in new coordinate system where Darcy velocity is parallel to x axis
p_rot  = jl.rotation_z(p,-θ) 
tp_rot = jl.rotation_z(tp,-θ) 

d = jl.evaluate_relevant_distances(jl.GroundWaterFlow(), p_rot, tp_rot) 
with np.nditer(d, op_flags=['readwrite']) as it:
    for x in it:
      y = x.item()
      if y[0] == 0.0 and y[1] == 0.0:
          x[...] = (0., params.rb, y[2], y[3])

tstep = 8760*3600/12.
tmax = 8760*3600*10.
t = jl.range(tstep, step = tstep, stop = tmax)
Nt = jl.length(t) # number of time steps


# 6. THERMAL RESPONSES
# mutual response function between pairs of segments (adiabatic surface boundary condition)
g = np.ones(d.shape + (len(t),))
with np.nditer(g, flags=['multi_index'], op_flags=['readwrite']) as it:
    for x in it:
      i = it.multi_index
      coord = d[i[0], i[1]]
      tt = t[i[2]]
      x[...] = 1 / (2*np.pi*params.λs)*jl.mfls_adiabatic_surface(tt, α, coord[0], coord[1], coord[2], vt, h, coord[3], atol = 1e-9)

# Matrix containing response function for each pair of segments at time-step 1
G = g[:,:,0]


# 7. LOADING CONDITION: for this particular simulation we impose temperature as boundary condition
T0 = 10.                                                 # undisturbed temperature
Tfin_constraint = np.array([90. if i%12 in range(0, 6) else 55. for i in range(0,Nt)]) # input temperature
# Tfin_constraint = 90*ones(Nt)


# 8. MATRIX ASSEMBLY EXAMPLE
M_injection = np.zeros((3*Nb+Ns,3*Nb+Ns))      # matrix describing topology of injection problem
M_extraction = np.zeros((3*Nb+Ns,3*Nb+Ns))     # matrix describing topology of extraction problem
b = np.zeros(3*Nb+Ns)                       # given term 
X = np.zeros((Nt, 3*Nb+Ns))                    # vector of unknowns
qprime = np.zeros((Nt, Ns))                  # heat injection per meter
Δqbcurrentsum = np.zeros(Ns)               # net heat injection on a given segment (this variable is needed by the solver)


jl.build_matrix_b(M_injection, Nb, Ns,  
              k_in, k_out, k_b,  
              G, 
              mf, cpf, bh_map, h,          
              internal_to_external 
            )

jl.build_matrix_b(M_extraction, Nb, Ns,
            k_in, k_out, k_b, 
            G,
            mf, cpf, bh_map, h,             
            external_to_internal
          )

jl.build_giventerm_b(b, Nb, Ns, Tfin_constraint[0], T0, internal_to_external)


def solve_problem_b(X, M_injection,M_extraction,
                    b,
                    Nb, Ns, Nt,
                    Tfin_constraint,
                    internal_to_external,external_to_internal,
                    qprime, g, T0,
                    Δqbcurrentsum, h):
    for i in range(0, Nt):
        M = M_injection if i%12 in range(0, 6) else M_extraction
        branches = internal_to_external if i%12 in range(0, 6) else external_to_internal

        jl.solve_full_convolution_step_b(X, M, b, i+1, Nb, Ns,
                        Tfin_constraint, branches,
                        qprime, g, T0,
                        Δqbcurrentsum, h                  
                    )
    

solve_problem_b(X,M_injection,M_extraction,
                b,
                Nb, Ns, Nt,
                Tfin_constraint,
                internal_to_external, external_to_internal,
                qprime, g, T0,
                Δqbcurrentsum, h                  
                )

# 9. EXTRACT FROM SOLUTION VECTOR X
Tfin  = X[:, 0:2*Nb-1:2]
Tfout = X[:, 1:2*Nb+1:2]
Tb    = X[:, 2*Nb:3*Nb] 
q     =  jl.cumsum(qprime, dims = 1) 

# output temperature to compute  energy and exergy echanged 
last_borehole_in_branch = np.array([[x[0][-1] if i%12 in range(0, 6) else x[1][-1] for i in range(0, Nt)] for x in zip(internal_to_external, external_to_internal)])
Tfos  = np.array([[Tfout[idx,i-1] for (idx,i) in enumerate(ll)] for ll in last_borehole_in_branch]).transpose()
Tfo   = np.mean(Tfos, axis=1)


# 10. KPIs COMPUTATION
nbranches = 8
Energy_exchanged   = mf * cpf * nbranches * (Tfin_constraint - Tfo) * (8760/12.)/1e6
Energy_exchanged_by_season = np.reshape(Energy_exchanged, (Nt//6, 6))
Energy_exchanged_by_season2 = np.sum(Energy_exchanged_by_season, axis=1)

E_injected   =  Energy_exchanged_by_season2[0::2]
E_extracted  =  Energy_exchanged_by_season2[1::2]

eta = - E_extracted / E_injected

Exergy_exchanged = mf * cpf * nbranches * (Tfin_constraint - Tfo - (T0 + 273.15) * np.log( (Tfin_constraint + 273.15) / (Tfo + 273.15) )  ) * (8760/12.)/1e6
Exergy_exchanged_by_season = np.reshape(Exergy_exchanged, (Nt//6, 6))
Exergy_exchanged_by_season2 = np.sum(Exergy_exchanged_by_season, axis=1)

Ex_injected   =  Exergy_exchanged_by_season2[0::2]
Ex_extracted  =  Exergy_exchanged_by_season2[1::2]

psi = - Ex_extracted / Ex_injected
