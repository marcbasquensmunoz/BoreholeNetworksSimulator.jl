using BoreholeNetworksSimulator
using CSV
using Statistics

function load_positions_from_file(file)
    data = CSV.read(file, values, header=true, decimal=',')
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
    network,            # Heat injection
    reverse(network)    # Heat extraction
]

Tf_injection = 40.
Tf_extraction = 5.

borehole_positions = load_positions_from_file(borehole_locations)

method = ConvolutionMethod()
medium = FlowInPorousMedium(λw = 0.6, λs = 2., Cw = 4.18*1e6, Cs = 1.7*1e6, θ = 0., Φ = 0.2, T0 = 10.)
borehole = SingleUPipeBorehole(H = 50., D = 4., λg = 2.5, pipe_position = ((0.03, 0.0), (-0.03, 0.0)))
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=borehole_positions)
constraint = uniform_InletTempConstraint([i%12 in 1:6 ? Tf_injection : Tf_extraction for i=1:Nt], n_branches(network))
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

function operator(i, Tin, Tout, Tb, q, configurations, mass_flow_containers)
    mf = 0.5
    active_network = configurations[i%12 in 1:6 ? 1 : 2]
    Nbr = n_branches(active_network)
    mass_flow_containers .= mf
    BoreholeOperation(active_network, @view mass_flow_containers[1:Nbr])
end

containers = @time initialize(options)
@time simulate!(operator=operator, options=options, containers=containers)

containers.X


############
# Draw plots
include("plots/borefield.jl")
include("plots/monitor_branch.jl")

monitored_branches = [3, 8]
color_ranges = [Pair(colorant"blue", colorant"darkorange"), Pair(colorant"red", colorant"green")]

plot_borefield(network, borehole_positions, distinguished_branches = monitored_branches, colors = color_ranges)
monitor_branch(containers, network.branches[monitored_branches[1]], color_ranges[1], options.t)


nbranch1, nbranch2 = 8, 3 
branch1 = network.branches[nbranch1]
branch2 = reverse(network).branches[nbranch2]
color_branch1 = range(colorant"blue", stop=colorant"darkorange", length=length(branch1));
color_branch2 = range(colorant"red", stop=colorant"green", length=length(branch2));

Nb = BoreholeNetworksSimulator.n_boreholes(network)
Tfin = containers.X[1:2:2*Nb, :]
Tfout = containers.X[2:2:2*Nb, :]
Tb = containers.X[2*Nb+1:3Nb, :]
Tfos  = hcat([[Tfout[idx,i] for (idx,i) in enumerate(ll)] for ll in BoreholeNetworksSimulator.first_boreholes(network)]...)
Tfo   = mean(Tfos, dims =2)
Tfo   = reshape(Tfo, length(Tfo))
q = containers.X[3*Nb+1:end, :]

# FIGURE 2: temperature and heat flow along selected boreholes

scene12 = Figure(size = (600, 450))
ax2 = scene12[1, 1] = Axis(scene12, ylabel = "T [°C]") #, xlabel = "time [months]")
for b in zip(branch1, color_branch1)
    @show b
    lines!(ax2, collect(options.t ./ 12Δt), Tfin[b[1], :], color = b[2], linewidth = 2.)
    # scatter!(ax2, collect(t ./ tstep), Tfin[:,b[1]], color = b[2], markersize = 12., marker = '▲')
end
#scatter!(ax2, collect(options.t ./ 12Δt), Tfos[nbranch1, :],  color = :black, markersize = 8.)
#scatter!(ax2, collect(options.t ./ 12Δt), Tfo,  color = :grey, markersize = 12.)

ax3 = scene12[2, 1] = Axis(scene12, ylabel = "q [W/m]", xlabel = "time [years]")
for b in zip(branch1, color_branch1)    
    lines!(ax3, collect(options.t ./ 12Δt), q[b[1], :], color = b[2], linewidth = 2.)
end

hidexdecorations!(ax2, grid = false)

scene12
