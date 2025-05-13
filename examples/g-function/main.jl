using BoreholeNetworksSimulator
using CairoMakie
include("make_figure.jl")
include("pygfunction_heat.jl")
include("pygfunction_mixed.jl")

H = 150.
D = 4.
rb = 0.075

mass_flow_per_branch = 0.25
α = 1e-6
k_s = 2.
k_g = 1.
k_p = 0.42
r_in = 0.015
r_out = 0.02
R_fp = 0.109

ts = H^2/(9α)           

pos = -0.05
pos1 = (pos, 0.)
pos2 = (0., pos)

spacings = [7.5, 15., 22.5, 30., 45., 1e15]
scenarios = [(2, 2), (4, 4)]

fig_unif_heat = make_figure(scenarios, spacings, compute_uniform_heat, "Uniform heat exchange rate", ts=ts, H=H, D=D, rb=rb)
fig_total_heat = make_figure(scenarios, spacings, compute_constant_total_heat, "Constant total heat exchange rate", ts=ts, H=H, D=D, rb=rb)

save("$(@__DIR__)/plots/uniform_heat.png", fig_unif_heat)
save("$(@__DIR__)/plots/mixed.png", fig_total_heat)