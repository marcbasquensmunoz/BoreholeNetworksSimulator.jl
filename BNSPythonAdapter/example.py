import src.adapter
from juliacall import Main as jl

Δt = 8760*3600/12.
Nt = 10*12

α = 1e-6
λ = 3.
T0 = 10.
medium = jl.GroundMedium(α=α, λ=λ, T0=T0)

D = 10.
H = 100.

borehole = jl.SingleUPipeBorehole(H=H, D=D)

σ = 5.
positions = jl.Array[jl.Tuple[jl.Float64, jl.Float64]]([(0., 0.), (0., σ)])
borefield = jl.EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)

configurations = [jl.BoreholeNetwork([[1], [2]])]

q1 = 5.
q2 = 5.
loads = jl.Array[jl.Float64]([q1, q2])
constraint = jl.constant_HeatLoadConstraint(loads, Nt)

options = jl.SimulationOptions(
    method = jl.ConvolutionMethod(),
    constraint = constraint,
    borefield = borefield,
    medium = medium,
    Δt = Δt,
    Nt = Nt,
    configurations = configurations
)

def operator(i, Tin, Tout, Tb, q, configurations):
    return jl.BoreholeOperation(configurations[0], jl.Array[jl.Float64]([2., 2.]))

containers = jl.initialize(options)
jl.simulate_b(operator=operator, options=options, containers=containers)

containers.X
