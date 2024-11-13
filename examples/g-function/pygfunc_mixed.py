import pygfunction as gt
import matplotlib.pyplot as plt
import numpy as np
import math

H = 150.
D = 4.
r_b = 0.075


mass_flow_per_branch = 0.25
alpha = 1.0e-6
k_s = 2.
k_g = 1.
k_p = 0.42
r_in = 0.015
r_out = 0.02
R_fp = gt.pipes.conduction_thermal_resistance_circular_pipe(r_in, r_out, k_p)

pos = -0.05
pos1 = (pos, 0.)
pos2 = (0., pos)

Nt = 100
dt = 3600*8760*100.

l = 7.5

#positions = [(0.,0.), (l,0.), (5*l,5*l)]
positions = [(i*l, j*l) for i in range(10) for j in range(10)]

n = len(positions)

boreholes = [gt.boreholes.Borehole(H=H, D=D, r_b=r_b, x=x[0], y=x[1]) for x in positions]
Utubes = [gt.pipes.SingleUTube(pos=[pos1, pos2], r_in=r_in, r_out=r_out,
                                borehole=bh, k_s=k_s, k_g=k_g, R_fp=R_fp, J=0) for bh in boreholes]

bore_connectivity = list(-1 * np.ones(len(positions)))
time = np.array([dt*(i+1) for i in range(Nt)])
ts = H**2/(9.*alpha)           
tts = np.log(time / ts)
cp_f = 4182.
network = gt.networks.Network(boreholes, Utubes, bore_connectivity, n * mass_flow_per_branch, cp_f)

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

T0 = 10.
Q = H * n   # Total heat extraction in the borefield per length equal to 1  

network = jl.all_parallel_network(n)

borehole = jl.SingleUPipeBorehole(H=H, D=D, rb=r_b, rpi=r_in, rpo=r_out, λg=k_g, pipe_position=(pos1, pos2))
borefield = jl.EqualBoreholesBorefield(borehole_prototype=borehole, positions=jl.Array[jl.Tuple[jl.Float64, jl.Float64]](positions))
medium = jl.GroundMedium(α=alpha, λ=k_s, T0=T0)
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

operator = jl.ConstantOperator(network, mass_flows = jl.Array[jl.Float64](mass_flow_per_branch * np.ones(n)))

containers = jl.initialize(options)
jl.simulate_b(operator=operator, options=options, containers=containers)

Tbm = (np.mean(containers.X[2*n:3*n, :], axis=0) - T0) * 2 * math.pi * k_s 

plt.figure(0)
plt.plot(tts, gfunc.gFunc)
plt.plot(tts, Tbm)

err = np.log10(abs(gfunc.gFunc - Tbm))

plt.figure(1)
plt.plot(tts, err)
