using BoreholeNetworksSimulator
using WGLMakie

function make_plot(axis, d)
    Δt = 3600.
    Nt = 8760*100

    D = 10.
    H = 10.

    n = 2
    m = 2
    Nb = n*m

    α = 1e-6
    λ = 3.

    Q = H

    network = all_parallel_network(n*m)
    configurations = [network]

    method = NonHistoryMethod()
    medium = GroundMedium(λ = λ, α = α, T0 = 0.)
    borehole = SingleUPipeBorehole(H = H, D = D, λg = 2.5, pipe_position = ((0.03, 0.0), (-0.03, 0.0)))
    borefield = RectangularBorefield(n, m, d, d, borehole)
    constraint = constant_HeatLoadConstraint(Q .* ones(n*m), Nt)
    fluid = Water()

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

    operator = SimpleOperator(mass_flow = 1., branches = n_branches(network))
    containers = @time initialize(options)

    @time simulate!(operator=operator, options=options, containers=containers)

    gfunc = sum(containers.X[2Nb+1:3Nb, :], dims=1) / Nb * (2π*λ)
    lines!(axis, log.(options.t ./ (H^2 / 9α)) , gfunc[1,:], label=L"\frac{d}{L} = %$(d/H)")
end


fig = Figure()
axis = fig[1, 1] = Axis(fig, ylabel = "g-function", xlabel = L"\ln \ \frac{t}{t_s}")

make_plot(axis, 0.05*10)
make_plot(axis, 0.1*10)
make_plot(axis, 0.25*10)
make_plot(axis, 0.5*10)

fig[1, 2] = Legend(fig, axis, "", framevisible = false)

fig
# save("examples/g-function/gfunction.png", fig)
