
import pygfunction as gt
import math 

dt = 3600.
alpha = 1e-6
kg = 3.

t = 3600*8760.

H1 = 100.
H2 = 150.
D1 = 50.
D2 = 0.

rb1 = 0.1
rb2 = rb1

sigma = 10.

borehole1 = gt.boreholes.Borehole(H=H1, D=D1, r_b=rb1, x=0., y=0.)
borehole2 = gt.boreholes.Borehole(H=H2, D=D2, r_b=rb2, x=0., y=sigma)

fls = gt.heat_transfer.finite_line_source(t, alpha, borehole1, borehole2) / (2 * math.pi * kg)



# -------------------------------------------------------------------------
# BoreholeNetworksSimulator version
# -------------------------------------------------------------------------

import sys, os
sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[0:-3]))
import BNSPythonAdapter.src.adapter
from juliacall import Main as jl
jl.seval("using FiniteLineSource")

atol = 1e-8
rtol = 0.

s = jl.FiniteLineSource.SegmentToSegment(D1=D1, D2=D2, H1=H1, H2=H2, σ=sigma)
s_image = jl.BoreholeNetworksSimulator.image(s)
params = jl.FiniteLineSource.Constants(Δt=dt, α=alpha, rb=rb1, kg=kg, b=10.)
fls2_real = jl.BoreholeNetworksSimulator.sts_response(s, params, t, atol=atol, rtol=rtol)
fls2_image = jl.BoreholeNetworksSimulator.sts_response(s_image, params, t, atol=atol, rtol=rtol)
fls2 = fls2_real - fls2_image

err = fls-fls2

print(err)