using BoreholeNetworksSimulator
using BNSPlots
using CSV
using Colors
using Parameters
using WGLMakie
using Statistics

function load_positions_from_file(file)
    data = CSV.read(file, values, header=true, decimal=',', delim=';')
    data = reduce(hcat, data)
    [(data[i, 2], data[i, 3]) for i in 1:size(data)[1]]
end

borehole_locations = "$(@__DIR__)/data/Braedstrup_borehole_coordinates.txt"
Δt = 3600.
Nt = 8760*10

Nb = 48
network = BoreholeNetwork(Nb)
connect_to_source!(network, [35, 43, 47, 44, 18, 19, 20, 27])
connect_in_series!(network, [35, 36, 29, 37, 30, 22])
connect_in_series!(network, [43, 48, 42, 41, 40, 34])
connect_in_series!(network, [47, 46, 45, 39, 32, 33])
connect_in_series!(network, [44, 38, 31, 24, 25, 26])
connect_in_series!(network, [18, 12, 11, 17, 16, 23])
connect_in_series!(network, [19, 13, 7, 6, 5, 10])
connect_in_series!(network, [20, 14, 8, 3, 2, 1])
connect_in_series!(network, [27, 28, 21, 15, 9, 4])
connect_to_sink!(network, [22, 34, 33, 26, 23, 10, 1, 4])

configurations = [  
    network,            # Heat extraction
    reverse(network)    # Heat injection
]

Tf_injection = 90.
Tf_extraction = 55.

borehole_positions = load_positions_from_file(borehole_locations)

year_in_dt = round(Int, 8670*3600 / Δt)
half_year_in_dt = round(Int, 4380*3600 / Δt)

method = OriginalNonHistoryMethod()
medium = GroundMedium(α=3 / (1.7*1e6), λ = 3., T0 = 10.)
borehole = SingleUPipeBorehole(H = 75., D = 4., λg = 1.5, pipe_position = ((0.03, 0.0), (-0.03, 0.0)))
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=borehole_positions)
constraint = uniform_InletTempConstraint([i%year_in_dt in 1:half_year_in_dt ? Tf_injection : Tf_extraction for i=1:Nt], n_branches(network))
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

@with_kw mutable struct IntermittentSeasonalOperator <: Operator
    mass_flows
    mass_flows_cont = zeros(length(mass_flows))
    seasonal_configuration
    active_steps
    inactive_steps
    currently_running::Bool = true # 0: not running, 1: running
    current_step = 0
end

function BoreholeNetworksSimulator.operate(operator::SeasonalOperator, i, options, X)
    active_network = options.configurations[operator.seasonal_configuration[i]]
    BoreholeOperation(network=active_network, mass_flows=operator.mass_flows)
end

function BoreholeNetworksSimulator.operate(operator::IntermittentSeasonalOperator, i, options, X)
    @unpack mass_flows, mass_flows_cont, seasonal_configuration, active_steps, inactive_steps, currently_running, current_step = operator
    active_network = options.configurations[seasonal_configuration[i]]

    if current_step >= (currently_running ? active_steps : inactive_steps) 
        currently_running = !currently_running
        current_step = 0
    end
    current_step += 1
    mass_flows_cont .= currently_running ? mass_flows : 0.

    @pack! operator = currently_running, current_step
    BoreholeOperation(network=active_network, mass_flows=mass_flows_cont)
end

seasonal_configuration = [i%year_in_dt in 1:half_year_in_dt ? 1 : 2 for i in 1:Nt]
constant_operator = SeasonalOperator(mass_flows=0.5 .* ones(n_branches(network)), seasonal_configuration=seasonal_configuration)
intermitent_operator = IntermittentSeasonalOperator(mass_flows=0.5 .* ones(n_branches(network)), seasonal_configuration=seasonal_configuration, active_steps = 12, inactive_steps = 12)
containers = @time initialize(options)

reset!(options)
@time simulate!(operator=constant_operator, options=options, containers=containers)
constant_sol = copy(containers)
reset!(options)
@time simulate!(operator=intermitent_operator, options=options, containers=containers)
intermitent_sol = copy(containers)

#########################################
# Borefield plots
#########################################

parallel_network = all_parallel_network(Nb)
parallel_network_plot = plot_borefield(parallel_network, borehole_positions)
save("$(@__DIR__)/plots/parallel_network.png", parallel_network_plot)

coloring = vcat([(i, colorant"red") for i in 1:16], [(i, colorant"blue") for i in 17:32], [(i, colorant"green") for i in 33:48])
groups_network_plot = plot_borefield(parallel_network, borehole_positions, distinguished_boreholes = coloring)
save("$(@__DIR__)/plots/groups_network.png", groups_network_plot)

series_network_plot = plot_borefield(network, borehole_positions)
save("$(@__DIR__)/plots/series_network.png", series_network_plot)

branch = boreholes_in_branch(network, first_bh=27)
colors = BNSPlots.make_color_range((colorant"navajowhite2", colorant"darkgreen"), length(branch)) 
branch_coloring = [(bh, colors[i]) for (i, bh) in enumerate(branch)]
legend_network_plot = plot_borefield(network, borehole_positions, distinguished_boreholes = branch_coloring, distinguished_marker_size = 24)
save("$(@__DIR__)/plots/legend_network.png", legend_network_plot)

#########################################
# Averaged plots
#########################################

# Averaging window
T = 24

average_int = zeros(4*Nb, Nt)
average_cont = zeros(4*Nb, Nt)
for t in 1:div(Nt, T)
    period = T*(t-1)+1:t*T
    average_int[:, period] .= mean(intermitent_sol.X[:, period], dims=2)
    average_cont[:, period] .= mean(constant_sol.X[:, period], dims=2)
end

intermitent_Tfin = BNSPlots.get_Tfin(containers) 
mf = intermitent_sol.mf[1, :]
flow = findall(>(0), mf)
average_Tfin = zeros(Nb, Nt)
for t in 1:div(size(average_Tfin, 2), T)
    period = T*(t-1)+1:t*T
    average_Tfin[:, period] .= mean(intermitent_Tfin[:, period[findall(in(flow), period)]], dims=2)
end

average_int[1:2:2Nb, :] .= average_Tfin

intermitent_sol_daily_average = SimulationContainers(X = average_int, M=zeros(4Nb, 4Nb), b=zeros(4Nb), mf=zeros(4Nb, 4Nb))
constant_sol_daily_average = SimulationContainers(X = average_cont, M=zeros(4Nb, 4Nb), b=zeros(4Nb), mf=zeros(4Nb, 4Nb))

average_intermitent = monitor(intermitent_sol_daily_average, boreholes_in_branch(network, first_bh=27), options.t, steps = 1:T:Nt, display = [:Tfin, :Tb, :q])
average_constant = monitor(constant_sol_daily_average, boreholes_in_branch(network, first_bh=27), options.t, steps = 1:T:Nt, display = [:Tfin, :Tb, :q])

max_q = maximum(BNSPlots.get_q(constant_sol_daily_average))
min_q = minimum(BNSPlots.get_q(constant_sol_daily_average))
ylims!(average_intermitent.layout.content[1].content.content[2].content, min_q, max_q)
ylims!(average_constant.layout.content[1].content.content[2].content, min_q, max_q)

max_T = 1.05 * max(maximum(BNSPlots.get_Tfin(constant_sol_daily_average)), maximum(BNSPlots.get_Tb(constant_sol_daily_average)))
min_T = 0.92 * min(minimum(BNSPlots.get_Tfin(constant_sol_daily_average)), minimum(BNSPlots.get_Tb(constant_sol_daily_average)))
ylims!(average_intermitent.layout.content[1].content.content[1].content, min_T, max_T)
ylims!(average_constant.layout.content[1].content.content[1].content, min_T, max_T)

save("$(@__DIR__)/plots/average_intermitent.png", average_intermitent)
save("$(@__DIR__)/plots/average_constant.png", average_constant)


total_heat_int = sum(BNSPlots.get_q(intermitent_sol_daily_average)[boreholes_in_branch(network, first_bh=27), 1:T:end], dims=1)*borehole.H
total_heat_const = sum(BNSPlots.get_q(constant_sol_daily_average)[boreholes_in_branch(network, first_bh=27), 1:T:end], dims=1)*borehole.H

first_bh = 27
last_bh = boreholes_in_branch(network, first_bh=first_bh)[end]
outlet_bh = [i%365 in 1:182 ? last_bh : first_bh for i in 1:div(Nt, 24)]#[i%year_in_dt in 1:half_year_in_dt ? last_bh : first_bh for i=1:Nt]
Tfout_int = BNSPlots.get_Tfout(intermitent_sol_daily_average)
Tfout_const = BNSPlots.get_Tfout(constant_sol_daily_average)
outlet_T_int = [Tfout_int[outlet_bh[div(i, T)+1], i] for i in 1:T:Nt]
outlet_T_const = [Tfout_const[outlet_bh[div(i, T)+1], i] for i in 1:T:Nt]

plot_range = 2900:3300
heat_comparison = Figure()
ax_T = Axis(heat_comparison[1, 1], title = "Outlet temperature", xlabel = "Time [days]", ylabel = "T [°C]")
lines!(ax_T, outlet_T_const[plot_range])
lines!(ax_T, outlet_T_int[plot_range])

ax_Q = Axis(heat_comparison[2, 1], title = "Total heat extraction", xlabel = "Time [days]", ylabel = "Total heat [kW]")
lines!(ax_Q, total_heat_const[plot_range] ./ 1000)
lines!(ax_Q, total_heat_int[plot_range] ./ 1000)

linkxaxes!(ax_T, ax_Q)
save("$(@__DIR__)/plots/total_heat_comparison.png", heat_comparison)

#########################################
# Zommed in plots
#########################################
get_period(year, week) = year*8760+(week-1)*24*7+1:year*8760+week*24*7
year = 9

observation_weeks = [1, 13, 26, 28, 39, 52]

for week in observation_weeks
    period = get_period(year, week)
    branch = boreholes_in_branch(network, first_bh=27)
    intermitent = monitor(intermitent_sol, branch, options.t, steps = period, display = [:Tfin, :Tb, :q])
    constant = monitor(constant_sol, branch, options.t, steps = period, display = [:Tfin, :Tb, :q])

    q_const = BNSPlots.get_q(constant_sol)[branch, period]
    q_int = BNSPlots.get_q(intermitent_sol)[branch, period]
    Tfin_const = BNSPlots.get_Tfin(constant_sol)[branch, period]
    Tfin_int = BNSPlots.get_Tfin(intermitent_sol)[branch, period]
    Tb_const = BNSPlots.get_Tb(constant_sol)[branch, period]
    Tb_int = BNSPlots.get_Tb(intermitent_sol)[branch, period]
    

    max_q = max(maximum(q_const), maximum(q_int))
    min_q = min(minimum(q_const), minimum(q_int))
    if max_q == 0. 
        max_q = -min_q * 0.1
    elseif min_q == 0.
        min_q = -max_q * 0.1
    end 
    if max_q > 0 
        max_q = max_q * 1.02
    else
        max_q = max_q * 0.95
    end
    if min_q > 0 
        min_q = min_q * 0.95
    else
        min_q = min_q * 1.02
    end

    max_T_const = max(maximum(Tfin_const), maximum(Tb_const))
    min_T_const = min(minimum(Tfin_const), minimum(Tb_const))
    max_T_int = max(maximum(Tfin_int), maximum(Tb_int))
    min_T_int = min(minimum(Tfin_int), minimum(Tb_int))
    max_T = max(max_T_const, max_T_int)
    min_T = min(min_T_const, min_T_int)
    if max_T > 0 
        max_T = max_T * 1.02
    else
        max_T = max_T * 0.95
    end
    if min_T > 0 
        min_T = min_T * 0.95
    else
        min_T = min_T * 1.02
    end

    ylims!(intermitent.layout.content[1].content.content[2].content, min_q, max_q)
    ylims!(constant.layout.content[1].content.content[2].content, min_q, max_q)

    ylims!(intermitent.layout.content[1].content.content[1].content, min_T, max_T)
    ylims!(constant.layout.content[1].content.content[1].content, min_T, max_T)

    save("$(@__DIR__)/plots/intermitent_year_$(year)_week_$(week).png", intermitent)
    save("$(@__DIR__)/plots/constant_year_$(year)_week_$(week).png", constant)
end
