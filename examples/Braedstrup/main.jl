using BoreholeNetworksSimulator
using BNSPlots
using CSV
using Colors
using Parameters
using WGLMakie

function load_positions_from_file(file)
    data = CSV.read(file, values, header=true, decimal=',', delim=';')
    data = reduce(hcat, data)
    [(data[i, 2], data[i, 3]) for i in 1:size(data)[1]]
end

borehole_locations = "$(@__DIR__)/data/Braedstrup_borehole_coordinates.txt"
Δt = 8760*3600/12.
Nt = 120

network = BoreholeNetwork([
    [35, 36, 29, 37, 30, 22],
    [43, 48, 42, 41, 40, 34],
    [47, 46, 45, 39, 32, 33],
    [44, 38, 31, 24, 25, 26],
    [18, 12, 11, 17, 16, 23],
    [19, 13, 7, 6, 5, 10],
    [20, 14, 8, 3, 2, 1],
    [27, 28, 21, 15, 9, 4]
])
configurations = [  
    network,            # Heat extraction
    reverse(network)    # Heat injection
]

Tf_injection = 90.
Tf_extraction = 55.

borehole_positions = load_positions_from_file(borehole_locations)

method = ConvolutionMethod()
medium = FlowInPorousMedium(λw = 0.6, λs = 2., Cw = 4.18*1e6, Cs = 1.7*1e6, θ = 0., Φ = 0.2, T0 = 10.)
borehole = SingleUPipeBorehole(H = 50., D = 4., λg = 1.5, pipe_position = ((0.03, 0.0), (-0.03, 0.0)))
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=borehole_positions)
constraint = uniform_InletTempConstraint([i%12 in 1:6 ? Tf_injection : Tf_extraction for i=1:Nt], n_branches(network))
fluid = Water()

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

@with_kw struct SeasonalOperator <: Operator
    mass_flows
    seasonal_configuration
end

function BoreholeNetworksSimulator.operate(operator::SeasonalOperator, i, options, X)
    active_network = options.configurations[operator.seasonal_configuration[i]]
    BoreholeOperation(active_network, operator.mass_flows)
end

operator = SeasonalOperator(mass_flows=0.5 .* ones(n_branches(network)), seasonal_configuration=[i%12 in 1:6 ? 2 : 1 for i in 1:Nt])
containers = @time initialize(options)

@time simulate!(operator=operator, options=options, containers=containers)

containers.X


############
# Draw plots

monitored_branches = [3, 8]
color_ranges = [(colorant"darkorange", colorant"blue"), (colorant"red", colorant"green")]

borefiled_plot = plot_borefield(network, borehole_positions, distinguished_branches = monitored_branches, colors = color_ranges)
branch1 = monitor(containers, network.branches[monitored_branches[1]], options.t, display = [:Tfin], color_pair=color_ranges[1])

# save("examples/Braedstrup/plots/Braedstrup_borefield.png", borefiled_plot)
# save("examples/Braedstrup/plots/branch1.png", branch1)

