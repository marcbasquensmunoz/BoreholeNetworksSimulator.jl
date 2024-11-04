using PythonCall
using CondaPkg
using BoreholeNetworksSimulator
using Statistics
using WGLMakie
CondaPkg.add("pygfunction")

gt = pyimport("pygfunction")
np = pyimport("numpy")

H = 150.
D = 4.
rb = 0.075

m_flow_network = 0.25
alpha = 1.0e-6
k_s = 2.
k_g = 1.
r_in = 0.015
r_out = 0.02
R_fp = 0.109

pos = -0.05
pos1 = (pos, 0.)
pos2 = (0., pos)

B = 7.5
n = 3
m = 2



function compute_uniform_heat(Nt, Δt, positions)
    ########################################################
    # pygfunction version
    ########################################################

    boreholes = pylist([gt.boreholes.Borehole(H=H, D=D, r_b=rb, x=x[1], y=x[2]) for x in positions])
    Utubes = pylist([gt.pipes.SingleUTube(pos=[pos1, pos2], r_in=r_in, r_out=r_out,
                                    borehole=bh, k_s=k_s, k_g=k_g, R_fp=R_fp) for bh in boreholes])

    bore_connectivity = -1 .* ones(length(positions))
    times = np.array(Δt .* collect(1:Nt))
    cp_f = 4182.
    network = gt.networks.Network(boreholes, Utubes, bore_connectivity, m_flow_network, cp_f)

    method = "detailed"
    options = Dict( "nSegments" => 1 )

    gfunc = gt.gfunction.gFunction(
        network, alpha, time=times, boundary_condition="UHTR",
        options=options, method=method)

    gfunc_res = pyconvert(Vector, gfunc.gFunc)


    ########################################################
    # BoreholeNetworksSimulator version
    ########################################################

    T0 = 10.
    Q = H 
    n = length(positions)

    network = all_parallel_network(n)

    borehole = SingleUPipeBorehole(H=H, D=D, rb=rb, rp=r_in, dpw = r_out-r_in, λg=k_g, pipe_position=(pos1, pos2))
    borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)
    medium = GroundMedium(α = alpha, λ = k_s, T0 = T0)
    constraint = uniform_HeatLoadConstraint([Q for i in range(1, Nt+1)], n_branches(network))
    method = NonHistoryMethod()
    fluid = EthanolMix()

    options = SimulationOptions(
        method = method,
        constraint = constraint,
        borefield = borefield,
        fluid = fluid,
        medium = medium,
        Δt = Δt,
        Nt = Nt,
        configurations = [network]
    )

    operator = ConstantOperator(network, mass_flows = m_flow_network .* ones(n))

    containers = initialize(options)
    simulate!(operator=operator, options=options, containers=containers)

    Tbm = (mean(containers.X[2*n+1:3*n, :], dims=1) .- T0)[:] * 2π * k_s 

    ts = H^2/(9alpha)           
    tts = @. log(options.t / ts)
    error_gfunc = abs.(gfunc_res - Tbm)

    return (tts, Tbm, error_gfunc)
end

########################################################
# Figure
########################################################

#scenarios = [(2, 2, B) for B in (7.5, 15., 22.5, 30., 45., 10000000.)]
scenarios = [(10, 10, B) for B in (7.5)]

fig = Figure()
grid = fig[1, 1] = GridLayout()

axis_gfunc = Axis(grid[1, 1], ylabel = L" g_{BNS}", title = "Uniform heat exchange rate")
axis_error = Axis(grid[2, 1], ylabel = L"\log_{10} \mid g_{pyg} - g_{BNS}\mid ", xlabel = L"\text{ln} \, \frac{t}{t_s}")

for (n, m, B) in scenarios
    positions = [(B*(i-1), B*(j-1)) for i in 1:n for j in 1:m]

    early = compute_uniform_heat(30, 3600*24., positions)
    mid = compute_uniform_heat(12*10, 3600*24*30., positions)
    late = compute_uniform_heat(100, 3600*8760*10., positions)
    far_late = compute_uniform_heat(20, 3600*8760*1000., positions)

    tts = vcat(early[1], mid[1], late[1], far_late[1])
    Tbm = vcat(early[2], mid[2], late[2], far_late[2])
    error_gfunc = vcat(early[3], mid[3], late[3], far_late[3])

    lines!(axis_gfunc, tts, Tbm)
    lines!(axis_error, tts, log10.(error_gfunc))
end

ylims!(axis_error, -15, 0)

linkxaxes!(axis_gfunc, axis_error)
hidexdecorations!(axis_gfunc, grid = false)
rowgap!(grid, 15)

fig

# save("examples/g-function/plots/uniform_heat.png", fig)