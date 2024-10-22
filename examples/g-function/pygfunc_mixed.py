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

positions = [(0.,0.), (l,0.), (5*l,5*l)]
#positions = [(i*l, j*l) for i in range(2) for j in range(2)]

boreholes = [gt.boreholes.Borehole(H=H, D=D, r_b=r_b, x=x[0], y=x[1]) for x in positions]
Utubes = [gt.pipes.SingleUTube(pos=[pos1, pos2], r_in=r_in, r_out=r_out,
                                borehole=bh, k_s=k_s, k_g=k_g, R_fp=R_fp) for bh in boreholes]

bore_connectivity = list(-1 * np.ones(len(positions)))
time = np.array([dt*(i+1) for i in range(Nt)])
ts = H**2/(9.*alpha)           
tts = np.log(time / ts)
cp_f = 4182.
network = gt.networks.Network(boreholes, Utubes, bore_connectivity, m_flow_network, cp_f)

method = 'detailed'
options = {'nSegments': 1}

gfunc = gt.gfunction.gFunction(
    network, alpha, time=time, boundary_condition='MIFT',
    options=options, method=method)


# -------------------------------------------------------------------------
# BoreholeNetworksSimulator version
# -------------------------------------------------------------------------

import sys, os
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-3]))
import BNSPythonAdapter.src.adapter
from juliacall import Main as jl
import numpy as np

n = len(positions)
T0 = 10.
Q = H * n   # Total heat extraction in the borefield per length equal to 1  

network = jl.all_parallel_network(n)

borehole = jl.SingleUPipeBorehole(H=H, D=D, rb=r_b, rp=r_in, dpw = r_out-r_in, λg=k_g, pipe_position=(pos1, pos2))
borefield = jl.EqualBoreholesBorefield(borehole_prototype=borehole, positions=jl.Array[jl.Tuple[jl.Float64, jl.Float64]](positions))
medium = jl.GroundMedium(α = alpha, λ = k_s, T0 = T0)
constraint = jl.TotalHeatLoadConstraint(jl.Array[jl.Float64]([Q for i in range(1, Nt+1)]))
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


# ------------------

a_in = -1.71574315e-001
a_b = -a_in

containers.M[0,:]

ki = containers.M[0,0]
ko = containers.M[0,1]
kb = containers.M[0,6]
k = ki+ko

kb_pred = - a_b * n*H / (cp_f * m_flow_network)
ki_pred = - a_in * n*H / (cp_f * m_flow_network) -1 
