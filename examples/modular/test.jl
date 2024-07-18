using BoreholeNetworksSimulator
using FiniteLineSource

Δt = 8760*3600/12.
Nt = 10
network = [
    BoreholeNetwork([[1]])
]
Nbr = length(network)

α = 1e-6
λ = 3.
rb = 0.1

positions = [(0., i) for i in 0:0]

T0 = 0.
borehole = SingleUPipeBorehole(H=1., D=0., rb=rb)
medium = GroundMedium(α=α, λ=λ)
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions, medium=medium, T0 = T0)
constraint = InletTempConstraint(10*ones(1))
#constraint = HeatLoadConstraint([1.])
fluid = Fluid(cpf = 4182., name = "INCOMP::MEA-20%")
function operator(i, Tin, Tout, Tb, q)
    BoreholeOperation(network[1], 0.1 .* ones(1))
end


options_c = SimulationOptions(
    method = ConvolutionMethod(),
    constraint = constraint,
    borefield = borefield,
    fluid = fluid,
    Δt = Δt,
    Nt = Nt
)

containers_c = initialize(options_c)
@time simulate(operator=operator, options=options_c, containers=containers_c)
X1 = containers_c.X


options_nh = SimulationOptions(
    method = NonHistoryMethod(),
    constraint = constraint,
    borefield = borefield,
    fluid = fluid,
    Δt = Δt,
    Nt = Nt
)
containers_nh = initialize(options_nh)
@time simulate(operator=operator, options=options_nh, containers=containers_nh)
X2 = containers_nh.X

X1-X2



### FLS
q = [1 for i=1:10]
I = zeros(length(q))

setup = FiniteLineSource.SegmentToSegment(D1=0., H1=1., D2=0., H2=1., σ=.1)
params = FiniteLineSource.Constants(Δt=Δt, rb=rb, b=2.)
precomp = FiniteLineSource.precompute_parameters(setup, params=params)
compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

BoreholeNetworksSimulator.sts(setup, params, tstep)
