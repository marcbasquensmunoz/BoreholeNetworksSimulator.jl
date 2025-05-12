using PythonCall
using CondaPkg
using BoreholeNetworksSimulator
using Statistics

CondaPkg.add("pygfunction")

gt = pyimport("pygfunction")
np = pyimport("numpy")

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
    method = ConvolutionMethod()
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
