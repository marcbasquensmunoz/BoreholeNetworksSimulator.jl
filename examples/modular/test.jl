using BoreholeNetworksSimulator
using FiniteLineSource

Δt = 8760*3600/12.
Nt = 1
configurations = [
    BoreholeNetwork([[1]])
]

α = 1e-6
λ = 3.
rb = 0.1

D = 0.
H = 1.

positions = [(0., 0.) for i in 0:0]

T0 = 0.
borehole = SingleUPipeBorehole(H=H, D=D, rb=rb)
medium = GroundMedium(α=α, λ=λ)
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions, T0 = T0)
#constraint = InletTempConstraint(10*ones(1))
constraint = HeatLoadConstraint([1.])
fluid = Fluid(cpf = 4182., name = "INCOMP::MEA-20%")
function operator(i, Tin, Tout, Tb, q)
    BoreholeOperation(configurations[1], 0.1 .* ones(1))
end


options_c = SimulationOptions(
    method = ConvolutionMethod(),
    constraint = constraint,
    borefield = borefield,
    fluid = fluid,
    medium = medium,
   #boundary_condition=NoBoundary(),
    Δt = Δt,
    Nt = Nt
)

containers_c = initialize(options_c)
@time simulate!(operator=operator, options=options_c, containers=containers_c)
X1 = containers_c.X


options_nh = SimulationOptions(
    method = NonHistoryMethod(),
    constraint = constraint,
    borefield = borefield,
    fluid = fluid,
    #boundary_condition=NoBoundary(),
    Δt = Δt,
    Nt = Nt
)
containers_nh = initialize(options_nh)
@time simulate!(operator=operator, options=options_nh, containers=containers_nh)
X2 = containers_nh.X

X1-X2



### FLS
q = [1 for i=1:Nt]
I = zeros(length(q))
params = FiniteLineSource.Constants(Δt=Δt, rb=rb, b=2.)

setup = FiniteLineSource.SegmentToSegment(D1=D, H1=H, D2=D, H2=H, σ=rb)
precomp = FiniteLineSource.precompute_parameters(setup, params=params)
compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

I2 = zeros(length(q))
setup2 = FiniteLineSource.SegmentToSegment(D1=-D, H1=-H, D2=D, H2=H, σ=rb)
precomp2 = FiniteLineSource.precompute_parameters(setup2, params=params)
compute_integral_throught_history!(setup2, I=I2, q=q, precomp=precomp2, params=params)

I-I2

