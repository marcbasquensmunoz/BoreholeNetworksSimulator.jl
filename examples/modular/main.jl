using GeometryTypes: Point2
using CSV, DataFrames, JLD2
using BTESGroundWaterSimulator

function load_borefield_from_file(file)
    df = CSV.File(file; decimal=',', delim = ';') |> DataFrame
    borehole_positions = [Point2(x, y) for (x,y) in zip(df.X,df.Y) ]
    EqualBoreholesBorefield(borehole_prototype=SingleUPipeBorehole(H=50., D=4.), positions=borehole_positions, medium=GroundWaterMedium())
end

function operator(i, Tin, Tout, Tb, Δq, Q)
    BoreholeOperation(networks[i%12 in 1:6 ? 2 : 1], 0.5 .* ones(8), 4182.)
end

networks = 
[
    [   
        [22,30,37,29,36,35], 
        [34,40,41,42,48,43],  
        [33,32,39,45,46,47],                        
        [26,25,24,31,38,44],  
        [23,16,17,11,12,18],  
        [10,5,6,7,13,19],     
        [1,2,3,8,14,20],      
        [4,9,15,21,28,27]                                         
    ]
    ,
    [
        [35,36,29,37,30,22],
        [43,48,42,41,40,34],
        [47,46,45,39,32,33],
        [44,38,31,24,25,26],
        [18,12,11,17,16,23],
        [19,13,7,6,5,10],
        [20,14,8,3,2,1],
        [27,28,21,15,9,4]
    ]
]
tstep = 8760*3600/12.
tmax  = 8760*3600*10.
Nt = div(tmax, tstep)

constraint = InletTempConstraint([i%12 in 1:6 ? 90. : 55. for i = 1:Nt])
cdir = @__DIR__
borehole_positions_file = "$cdir/../example1/data/Braedstrup_borehole_coordinates.txt"
borefield = load_borefield_from_file(borehole_positions_file)

cache = ""
cache = "$cdir/results/simulation/cache_3.1536e7.jld2"
cache = "$cdir/results/simulation/cache_3.1536e8.jld2"


parameters = compute_parameters(borefield=borefield, tstep=tstep, tmax=tmax)
model = ConvolutionGroundModel(T0 = 10., parameters=parameters)
containers = SimulationContainers(parameters)
load_cache!(containers=containers, parameters=parameters, cache=cache)

@btime simulate(parameters=parameters, containers=containers, operator=operator, borefield=borefield, constraint=constraint, model=model)

save_cache(containers=containers, parameters=parameters, path=@__DIR__, title="test")
