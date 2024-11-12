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

mass_flow_per_branch = 0.25
α = 1.0e-6
k_s = 2.
k_g = 1.
r_in = 0.015
r_out = 0.02
R_fp = 0.109

ts = H^2/(9α)           

pos = -0.05
pos1 = (pos, 0.)
pos2 = (0., pos)


function compute_uniform_heat(Nt, Δt, positions)
    n = length(positions)

    ########################################################
    # pygfunction version
    ########################################################

    boreholes = pylist([gt.boreholes.Borehole(H=H, D=D, r_b=rb, x=x[1], y=x[2]) for x in positions])
    Utubes = pylist([gt.pipes.SingleUTube(pos=[pos1, pos2], r_in=r_in, r_out=r_out,
                                    borehole=bh, k_s=k_s, k_g=k_g, R_fp=R_fp) for bh in boreholes])

    bore_connectivity = -1 .* ones(length(positions))
    times = np.array(Δt .* collect(1:Nt))
    cp_f = 4182.
    network = gt.networks.Network(boreholes, Utubes, bore_connectivity, n * mass_flow_per_branch, cp_f)

    method = "detailed"
    options = Dict( "nSegments" => 1 )

    gfunc = gt.gfunction.gFunction(
        network, α, time=times, boundary_condition="UHTR",
        options=options, method=method)

    gfunc_res = pyconvert(Vector, gfunc.gFunc)


    ########################################################
    # BoreholeNetworksSimulator version
    ########################################################

    T0 = 10.
    Q = H 

    network = all_parallel_network(n)

    borehole = SingleUPipeBorehole(H=H, D=D, rb=rb, rpi=r_in, rpo=r_out, λg=k_g, pipe_position=(pos1, pos2))
    borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=positions)
    medium = GroundMedium(α = α, λ = k_s, T0 = T0)
    constraint = uniform_HeatLoadConstraint([Q for i in range(1, Nt+1)], n_branches(network))
    method = ConvolutionMethod()#NonHistoryMethod()
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

    operator = ConstantOperator(network, mass_flows = mass_flow_per_branch .* ones(n))

    containers = initialize(options)
    simulate!(operator=operator, options=options, containers=containers)

    Tbm = (mean(containers.X[2*n+1:3*n, :], dims=1) .- T0)[:] * 2π * k_s 

    tts = @. log(options.t / ts)
    error_gfunc = abs.(gfunc_res - Tbm)

    return (tts, Tbm, error_gfunc)
end

########################################################
# Figure
########################################################

scenarios = (2, 2, [7.5, 15., 22.5, 30., 45., 1e15])


n = scenarios[1]
m = scenarios[2]
BB = scenarios[3]

fig = Figure()
grid = fig[1, 1] = GridLayout()

axis_gfunc = Axis(grid[1, 1], ylabel = L" g_{BNS}", title = "Uniform heat exchange rate; $n x $m grid")
axis_error = Axis(grid[2, 1], ylabel = L"\log_{10} \mid g_{pyg} - g_{BNS}\mid ", xlabel = L"\text{ln} \, \frac{t}{t_s}")

for B in BB
    positions = [(B*(i-1), B*(j-1)) for i in 1:n for j in 1:m]

    early = compute_uniform_heat(30, 3600*24., positions)
    mid = compute_uniform_heat(12*10, 3600*24*30., positions)
    late = compute_uniform_heat(100, 3600*8760*10., positions)
    far_late = compute_uniform_heat(20, 3600*8760*1000., positions)

    tts = vcat(early[1], mid[1], late[1], far_late[1])
    Tbm = vcat(early[2], mid[2], late[2], far_late[2])
    error_gfunc = vcat(early[3], mid[3], late[3], far_late[3])

    zero_error = findall(x -> x==0., error_gfunc)
    error_gfunc[zero_error] .= eps()

    legend_value = B > 1e8 ? "\\infty" : B/H
    lines!(axis_gfunc, tts, Tbm, label = L"%$(legend_value)")
    lines!(axis_error, tts, log10.(error_gfunc))
end

ylims!(axis_error, -18, 0)

linkxaxes!(axis_gfunc, axis_error)
hidexdecorations!(axis_gfunc, grid = false)
rowgap!(grid, 15)

rbH = rb/H
DH = round(D/H, digits=2)

legend = fig[1, 2] = GridLayout()

Label(legend[1, 1], L"t_s = 2.5 \times 10^{9} \ \text{s}")
Label(legend[2, 1], L"\frac{r_b}{H} = %$(rbH)")
Label(legend[3, 1], L"\frac{D}{H} = %$(DH)")
Legend(legend[4, 1], axis_gfunc, L"B/H")

fig

save("examples/g-function/plots/uniform_heat_$(n)x$(m).png", fig)