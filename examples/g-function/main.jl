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

    network = BoreholeNetwork([[i] for i in 1:n*m])
    configurations = [network]

    function create_rectangular_field(n, m, d)
        [((i-1)*d, (j-1)*d) for i in 1:n for j in 1:m]
    end

    method = NonHistoryMethod()
    medium = GroundMedium(λ = λ, α = α, T0 = 20.)
    borehole = SingleUPipeBorehole(H = H, D = D, λg = 2.5, pipe_position = ((0.03, 0.0), (-0.03, 0.0)))
    borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=create_rectangular_field(n, m, d))
    constraint = constant_HeatLoadConstraint(ones(n*m), Nt)
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

    operator = SimpleOperator(mass_flow = 1., branches =  n_branches(network))
    containers = @time initialize(options)

    @time simulate!(operator=operator, options=options, containers=containers)

    bh = Int(2Nb+1 + Nb/2)
    lines!(axis, log.(options.t ./ (H^2 / 9α)) , containers.X[bh, :], label=L"\frac{d}{L} = %$(d/H)")
end


fig = Figure()
axis = fig[1, 1] = Axis(fig, ylabel = "g-function", xlabel = L"\ln \ \frac{t}{t_s}")

make_plot(axis, 0.05*10)
make_plot(axis, 0.1*10)
make_plot(axis, 0.25*10)
make_plot(axis, 0.5*10)

fig[1, 2] = Legend(fig, axis, "", framevisible = false)

fig
save("examples/g-function/gfunction.png", fig)
