import pygfunction as gt
import matplotlib.pyplot as plt
import numpy as np
import math

H = 150.
D = 4.
r_b = 0.075


m_flow_network = 0.25
alpha = 1.0e-6
k_s = 2.
k_g = 1.
r_in = 0.015
r_out = 0.02
R_fp = 0.109

pos = -0.05
pos1 = (pos, 0.)
pos2 = (0., pos)

Nt = 12*100
dt = 3600*24*30.

l = 10.

# High error: 
#positions = [(0.,0.), (l,0.), (5*l,5*l)]
positions = [(i*l, j*l) for i in range(2) for j in range(2)]

boreholes = [ gt.boreholes.Borehole(H=H, D=D, r_b=r_b, x=x[0], y=x[1]) for x in positions]
Utubes = [gt.pipes.SingleUTube(pos=[pos1, pos2], r_in=r_in, r_out=r_out,
                                borehole=bh, k_s=k_s, k_g=k_g, R_fp=R_fp) for bh in boreholes]

bore_connectivity = list(-1 * np.ones(len(positions)))
time = np.array([dt*(i+1) for i in range(Nt)])
ts = H**2/(9.*alpha)           
tts = np.log(time / ts)
cp_f = 4182.
network = gt.networks.Network(boreholes, Utubes, bore_connectivity, m_flow_network, cp_f)

method = 'detailed'#'similarities'
options = {'nSegments': 1}

gfunc = gt.gfunction.gFunction(
    network, alpha, time=time, boundary_condition='UHTR',
    options=options, method=method)


# -------------------------------------------------------------------------
# BoreholeNetworksSimulator version
# -------------------------------------------------------------------------

import sys, os
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-3]))
import BNSPythonAdapter.src.adapter
from juliacall import Main as jl
jl.seval("using FiniteLineSource")
import numpy as np
import pandas as pd

n = len(positions)
T0 = 10.
Q = H 

network = jl.all_parallel_network(n)

borehole = jl.SingleUPipeBorehole(H=H, D=D, rb=r_b, rp=r_in, dpw = r_out-r_in, λg=k_g, pipe_position=(pos1, pos2))
borefield = jl.EqualBoreholesBorefield(borehole_prototype=borehole, positions=jl.Array[jl.Tuple[jl.Float64, jl.Float64]](positions))
medium = jl.GroundMedium(α = alpha, λ = k_s, T0 = T0)
constraint = jl.uniform_HeatLoadConstraint(jl.Array[jl.Float64]([Q for i in range(1, Nt+1)]), jl.n_branches(network))
method = jl.ConvolutionMethod()
fluid = jl.EthanolMix()

options = jl.SimulationOptions(
    method = method,
    constraint = constraint,
    borefield = borefield,
    fluid = fluid,
    medium = medium,
    Δt = dt,
    Nt = Nt,
    configurations = jl.Array[jl.BoreholeNetwork]([network])
)

operator = jl.SimpleOperator(mass_flow = m_flow_network, branches = jl.n_branches(network))

containers = jl.initialize(options)
jl.simulate_b(operator=operator, options=options, containers=containers)

Tbm = (np.mean(containers.X[2*n:3*n, :], axis=0) - T0) * 2 * math.pi * k_s 

plt.figure(0)
plt.plot(tts, gfunc.gFunc)
plt.plot(tts, Tbm)

err = abs(gfunc.gFunc - Tbm)

plt.figure(1)
plt.plot(tts, err)


# --------------
# Heat transfer 
# --------------
#M_py = np.matrix([[gt.heat_transfer.finite_line_source(dt, alpha, boreholes[i], boreholes[j]) / (2 * math.pi * k_s) for i in range(4)] for j in range(4)])
#M_jl = containers.M[8:12,12:16]
#print(M_py - M_jl)


# ---------------
# Internal model
# ---------------
kk = Utubes[0].coefficients_outlet_temperature(m_flow_network, cp_f, 1)
# T_{f,out}} = a_{in} T_{f,in} + a_{b} T_b
#  pipes.py, ln 1066: def coefficients_outlet_temperature

kin_jl = containers.M[0,0]
kout_jl = containers.M[0,1]
kb_jl = containers.M[0,8]

print(kk[0] - (-kin_jl/kout_jl))
print(kk[1] - (-kb_jl/kout_jl))



# ---------------
# Try to compute g-func
# ---------------

def get_fls(t, D1, D2, H1, H2, sigma, rb, k_s, alpha, dt):
    atol = 0.
    rtol = 1e-8
    params = jl.FiniteLineSource.Constants(Δt=dt, α=alpha, rb=rb, kg=k_s, b=10.)
    s = jl.FiniteLineSource.SegmentToSegment(D1=D1, D2=D2, H1=H1, H2=H2, σ=sigma)
    fls = jl.BoreholeNetworksSimulator.response(jl.DirichletBoundaryCondition(), s, params, t, atol=atol, rtol=rtol)
    return fls

def get_dist(b1, b2):
    if b1.x == b2.x and b1.y == b2.y:
        return b1.r_b 
    else:
        return np.sqrt((b1.x-b2.x)**2+(b1.y-b2.y)**2)

def get_g_func(t):
    M = np.matrix([[get_fls(t, D, D, H, H, get_dist(boreholes[i], boreholes[j]), r_b, k_s, alpha, dt) for i in range(3)] for j in range(3)])
    Tb = np.sum(M, axis=1)
    return np.mean(Tb) * 2 * k_s * math.pi

g_func = [get_g_func(dt*(i+1)) for i in range(Nt)]

err2 = g_func - gfunc.gFunc