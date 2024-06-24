using BoreholeNetworksSimulator

tstep = 8760*3600/12.
tmax  = 8760*3600/12.
Nt = Int(div(tmax, tstep))

network = [[1, 2]]
borefield = EqualBoreholesBorefield(borehole_prototype=SingleUPipeBorehole(H=50., D=0.), positions=[(0., 0.), (1., 1.)], medium=GroundMedium(), T0 = 10.)
parameters = compute_parameters(borefield=borefield, tstep=tstep, tmax=tmax)
constraint = InletTempConstraint(90*ones(Nt))
method = NonHistoryMethod(parameters=parameters, borefield=borefield)
containers = SimulationContainers(parameters)
function operator(i, Tin, Tout, Tb, Δq, Q)
    BoreholeOperation(network, ones(1), 4182.)
end

simulate(parameters=parameters, containers=containers, operator=operator, borefield=borefield, constraint=constraint, method=method)