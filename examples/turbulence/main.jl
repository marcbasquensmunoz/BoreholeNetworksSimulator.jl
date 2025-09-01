using BoreholeNetworksSimulator
using BNSPlots
using BNSPlots: get_Tfin, get_Tfout, get_Tb, get_q, get_mf
using Graphs
using WGLMakie

Δt = 3600.
Nt = 8760

D = 0.
H = 100.

Nb = 1

α = 1e-6
λ = 3.

laminar_mf_low = 0.22
laminar_mf = 0.23
turbulent_mf = 0.24

network = all_parallel_network(Nb)
positions = [(0., 0.)]
configurations = [network]

method = OriginalNonHistoryMethod()
medium = GroundMedium(λ=λ, α=α, T0=9.)
borehole = SingleUPipeBorehole(H=H, D=D, pipe_position = ((0.03, 0.0), (-0.03, 0.0)))
borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)
constraint = constant_InletTempConstraint(10 .* ones(Nb), Nt)
fluid = EthanolMix()

options = SimulationOptions(
    method = method,
    constraint = constraint,
    borefield = borefield,
    fluid = fluid,
    medium = medium,
    boundary_condition = DirichletBoundaryCondition(),
    Δt = Δt,
    Nt = Nt,
    configurations = configurations
)

struct VariableMassFlowOperator <: Operator 
    mf1
    mf2
end

function BoreholeNetworksSimulator.operate(op::VariableMassFlowOperator, step, options, X)
    mass_flow = floor((step / (24*30)))%2 == 0 ? op.mf1 : op.mf2
    network = options.configurations[1]
    source_valve = absolute_valve(Graphs.outneighbors(network.graph, source(network)), [mass_flow])
    BoreholeOperation(Dict(source(network) => source_valve), mass_flow, network)
end

operator_phase_transition = VariableMassFlowOperator(laminar_mf, turbulent_mf)
operator_laminar = VariableMassFlowOperator(laminar_mf_low, laminar_mf)

containers = @time initialize(options)
@time simulate!(operator=operator_phase_transition, options=options, containers=containers)
sol_phase_transition = copy(containers)
reset!(options)
@time simulate!(operator=operator_laminar, options=options, containers=containers)
sol_laminar = copy(containers)

# Plot
steps = 1:Nt
ticks = 0:12
time_axis = steps ./ (24*30)
borehole = 1
color_laminar = Makie.wong_colors()[1]
color_transition = Makie.wong_colors()[2]
fig = Figure()
axis_T = Axis(fig[1, 1], ylabel = L"T_{f, \text{out}} \, \left[ °C \right]")
axis_Q = Axis(fig[2, 1], ylabel = L"q \, \left[ \frac{W}{m} \right]")
axis_mf = Axis(fig[3, 1], ylabel = L"\dot{m} \, \left[ \frac{kg}{s} \right]")
axis_mf.xlabel = "time [months]"
axis_T.xticks = ticks
axis_Q.xticks = ticks
axis_mf.xticks = ticks

lines!(axis_T, time_axis, get_Tfout(sol_laminar)[borehole, steps], color = color_laminar, linewidth = 2., label="Laminar flow")
lines!(axis_T, time_axis, get_Tfout(sol_phase_transition)[borehole, steps], color = color_transition, linewidth = 2., label="Phase transition")

lines!(axis_Q, time_axis, get_q(sol_laminar)[borehole, steps], color = color_laminar, linewidth = 2.)
lines!(axis_Q, time_axis, get_q(sol_phase_transition)[borehole, steps], color = color_transition, linewidth = 2.)

lines!(axis_mf, time_axis, get_mf(sol_laminar)[borehole, steps], color = color_laminar, linewidth = 2.)
lines!(axis_mf, time_axis, get_mf(sol_phase_transition)[borehole, steps], color = color_transition, linewidth = 2.)

hidexdecorations!(axis_T, grid = false)
hidexdecorations!(axis_Q, grid = false)

linkxaxes!(axis_T, axis_Q)
linkxaxes!(axis_Q, axis_mf)

axislegend(axis_T, merge = true, unique = true, position=:rb, orientation=:horizontal)
fig
save("$(@__DIR__)/plots/comparison.png", fig)

#fig_laminar = monitor(sol_laminar, [1], options.t, display = [:Tfout, :Tb, :q, :mf])
#fig_transition = monitor(sol_phase_transition, [1], options.t, display = [:Tfout, :Tb, :q, :mf])
