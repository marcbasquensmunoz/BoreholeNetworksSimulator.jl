using BoreholeNetworksSimulator
import BoreholeNetworksSimulator: n_branches
using CSV

function load_borefield_from_file(file)
    data = CSV.read(file, values, header=true)
    data = reduce(hcat, data)
    [(data[i, 2], data[i, 3]) for i in 1:size(data)[1]]
end

borehole_locations = "$(@__DIR__)/Braedstrup_borehole_coordinates.txt"
Δt = 8760*3600/12.
Nt = 120

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
configurations = [
    network,
    reverse(network)
]

method = ConvolutionMethod()
medium = FlowInPorousMedium(λw = 0.6, λs = 2., Cw = 4.18*1e6, Cs = 1.7*1e6, θ = 0., Φ = 0.2, T0 = 10.)
borehole = SingleUPipeBorehole(H = 50., D = 4., λg = 2.5, pipe_position = ((0.03, 0.0), (-0.03, 0.0)))
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=load_borefield_from_file(borehole_locations))
constraint = uniform_InletTempConstraint([i%12 in 1:6 ? 90. : 55. for i=1:Nt], n_branches(network))
fluid = Fluid(cpf = 4182., name = "INCOMP::MEA-20%")

options = SimulationOptions(
    method = method,
    constraint = constraint,
    borefield = borefield,
    fluid = fluid,
    medium = medium,
    Δt = Δt,
    Nt = Nt,
    configurations = configurations
)

function operator(i, Tin, Tout, Tb, q, configurations)
    mf = 0.5
    active_network = configurations[i%12 in 1:6 ? 1 : 2]
    Nbr = n_branches(active_network)
    BoreholeOperation(active_network, mf .* ones(Nbr))
end

containers = @time initialize(options)
@time simulate!(operator=operator, options=options, containers=containers)

containers.X