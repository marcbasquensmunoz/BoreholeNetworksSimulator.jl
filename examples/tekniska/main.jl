using BoreholeNetworksSimulator
using WGLMakie
using BNSPlots
using CSV
using Statistics

Δt = 3600.
Nt = 24*7

network = BoreholeNetwork([[i] for i in 1:10])
configurations = [  
    network
]


σ = 10.
Δy = 5.
borehole_positions = vcat([(σ*(i-1), 0.) for i in 1:5], [(σ*(i-1/2), Δy) for i in 1:5])
Q_ref = 40000.

affluency = [
    [0.2, 0.33, 0.4, 0.47, 0.47, 0.45, 0.4, 0.37, 0.4, 0.47, 0.45, 0.35],
    [0.25, 0.35, 0.5, 0.56, 0.55, 0.5, 0.35, 0.3, 0.25, 0.2, 0.1, 0.05],
    [0.28, 0.38, 0.55, 0.6, 0.63, 0.58, 0.42, 0.37, 0.33, 0.2, 0.12, 0.08],
    [0.33, 0.45, 0.63, 0.7, 0.7, 0.65, 0.57, 0.45, 0.33, 0.25, 0.15, 0.08],
    [0.28, 0.4, 0.55, 0.6, 0.6, 0.57, 0.57, 0.6, 0.66, 0.7, 0.66, 0.53],
    [0.33, 0.57, 0.75, 0.9, 0.85, 0.73, 0.66, 0.62, 0.6, 0.6, 0.6, 0.5],
    [0.25, 0.4, 0.6, 0.66, 0.65, 0.57, 0.5, 0.33, 0.25, 0.17, 0.1, 0.05]
]

function plot_affluency_day(i)
    fig = Figure()
    ax = Axis(fig[1, 1], yreversed = true)  
    scatter!(ax, affluency[i])
    ylims!(ax, (0, 1))
    fig
end

Q_week = [i%24 in 10:21 ? Q_ref * affluency[ceil(Int, i/(24))][i%24-9] : 0. for i in 1:24*7]
seasonal_shift = 0.6 * [sin(2π*i/52) + 1.5 for i in 1:52] 
Q_year = collect(Iterators.flatten([T .* Q_week for T in seasonal_shift]))


function plot_Q()
    fig = Figure(size=(1500, 400))
    ax = Axis(fig[1, 1])  
    lines!(ax, Q)
    fig
end

#method = ConvolutionMethod()
method = NonHistoryMethod()
medium = GroundMedium(λ = 3.1, α = 1e-6, T0 = 9.)
borehole = SingleUPipeBorehole(H = 240., D = 4., λg = 2.5, rb=0.145/2, pipe_position = ((0.03, 0.0), (-0.03, 0.0)))
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=borehole_positions)
constraint = TotalHeatLoadConstraint(Q_year)
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
    mf1 = 0.5
    mf2 = 0.5
    active_network = configurations[1]
    Nbr = n_branches(active_network)
    mass_flow_containers[1:6] .= mf1
    mass_flow_containers[7:10] .= mf2
    BoreholeOperation(active_network, @view mass_flow_containers[1:Nbr])
end

containers = @time initialize(options)
@time simulate!(operator=operator, options=options, containers=containers)


plot_borefield(network, borehole_positions)
monitor(containers, [1, 2], options.t)

Q_tot = sum(containers.X[31:40, :], dims=1)[:]
lines(240*Q_tot)
