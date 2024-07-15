using BoreholeNetworksSimulator
using FiniteLineSource

tstep = 8760*3600/12.
tmax  = 10*8760*3600/12.
Nt = Int(div(tmax, tstep))
network = [
    BoreholeNetwork([[1, 2], [3]])
]
Nbr = length(network)

α = 1e-6
λ = 3.
rb = 0.1

positions = [(0., i) for i in 0:2]


borefield = EqualBoreholesBorefield(borehole_prototype=SingleUPipeBorehole(H=10., D=0., rb=rb), positions=positions, medium=GroundMedium(α=α, λ=λ), T0 = 0.)
parameters = compute_parameters(borefield=borefield, tstep=tstep, tmax=tmax)
constraint = InletTempConstraint(10*ones(Nt))
#constraint = HeatLoadConstraint([1., 1.])
function operator(i, Tin, Tout, Tb, q)
    BoreholeOperation(network[1], 0.1 .* ones(Nbr), 4182.)
end


convolution = ConvolutionMethod(parameters=parameters, borefield=borefield)
containers_c = SimulationContainers(parameters)
@time simulate(parameters=parameters, containers=containers_c, operator=operator, borefield=borefield, constraint=constraint, method=convolution)
X1 = containers_c.X

nonhistory = NonHistoryMethod(parameters=parameters, borefield=borefield, b = 10.)
containers_nh = SimulationContainers(parameters)
@time simulate(parameters=parameters, containers=containers_nh, operator=operator, borefield=borefield, constraint=constraint, method=nonhistory)
X2 = containers_nh.X

X1-X2



### FLS
q = [1 for i=1:10]
I = zeros(length(q))
#q = X1[4,:]
Δt = tstep

setup = FiniteLineSource.SegmentToSegment(D1=0., H1=1., D2=0., H2=1., σ=1.)
params = FiniteLineSource.Constants(Δt=Δt, rb=rb, b=2.)
precomp = FiniteLineSource.precompute_parameters(setup, params=params)
compute_integral_throught_history!(setup, I=I, q=q, precomp=precomp, params=params)

BoreholeNetworksSimulator.sts(setup, params, tstep)


function send_result_to_database(i) sleep(0.1) end
function update_database_result_count(nresults) end

function test()
    simulation_results = collect(1:10)
    @sync for result in simulation_results
        Dagger.@spawn occupancy=Dict(Dagger.ThreadProc=>0.) send_result_to_database(result)
    end
    nresults = length(simulation_results)
    wait(Dagger.@spawn update_database_result_count(nresults))
    update_database_result_count(nresults)
end

