using GeometryTypes: Point2
using CSV, DataFrames, JLD2
using BoreholeNetworksSimulator

function load_borefield_from_file(file)
    df = CSV.File(file; decimal=',', delim = ';') |> DataFrame
    borehole_positions = [Point2(x, y) for (x,y) in zip(df.X,df.Y)]
    EqualBoreholesBorefield(borehole_prototype=SingleUPipeBorehole(H=50., D=4.), positions=borehole_positions, medium=GroundMedium(), T0 = 10.)
end

network = BoreholeNetwork([
    [22,30,37,29,36,35], 
    [34,40,41,42,48,43],  
    [33,32,39,45,46,47],                        
    [26,25,24,31,38,44],  
    [23,16,17,11,12,18],  
    [10,5,6,7,13,19],     
    [1,2,3,8,14,20],      
    [4,9,15,21,28,27] 
])

configurations = [network, reverse(network)]

tstep = 8760*3600/12.
tmax  = 8760*3600.
Nt = div(tmax, tstep)

cdir = @__DIR__
borehole_positions_file = "$cdir/../old/data/Braedstrup_borehole_coordinates.txt"
borefield = load_borefield_from_file(borehole_positions_file)
parameters = compute_parameters(borefield=borefield, tstep=tstep, tmax=tmax)
constraint = InletTempConstraint(90 .* ones(8))
method = ConvolutionMethod(parameters=parameters, borefield=borefield)
#method = NonHistoryMethod(parameters=parameters, borefield=borefield)
containers = SimulationContainers(parameters)

function operator(i, Tin, Tout, Tb, Q)
    BoreholeOperation(configurations[i%12 in 1:6 ? 2 : 1], 0.5 .* ones(8), 4182.)
end

simulate(parameters=parameters, containers=containers, operator=operator, borefield=borefield, constraint=constraint, method=method)
save_cache(containers=containers, parameters=parameters, path=@__DIR__, title="test")
