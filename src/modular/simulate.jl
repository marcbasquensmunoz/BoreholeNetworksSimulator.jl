# X contains: Tin     | Tout    | Tb        | q
#             1:2:2Nb | 2:2:2Nb | 2Nb+1:3Nb | 3Nb+1:3Nb+Ns
#-----------------------------------------------------------
# Amount:       Nb    |   Nb    |    Nb     |    Ns

########################################################
# Equations present in the system
#   (1)  Internal model
#   (2)  Branch topology and temperature constraints
#   (3)  Heat transfer model
#   (4)  Heat balance equations

"""
    simulate!(;options::SimulationOptions, operator, containers::SimulationContainers)

Run the simulation defined by `options`. 
At the end of simulation, `containers.X` will contain the results. `containers` should be the output of [`initialize`](@ref).

`operator` should be a function that returns a `BoreholeOperation` and with signature `operator(step, options, X)`:
- `step::Int` is the time step.
- `options` is the `SimulationOptions` object passed to `simulate!`.
- `X` is the matrix `containers.X` containing `Tin`, `Tout`, `Tb`, and `q` for each borehole.
"""
function simulate!(;operator, options::SimulationOptions, containers::SimulationContainers)
    @unpack medium = options
    @unpack Nb, Ts, Nt = options
    
    mass_flows = zeros(Nb+2)
    fluid_T = get_T0(medium) .* ones(2Nb)
    node_queue = Queue{Tuple{Int, Float64}}()

    # Simulation loop
    for i = Ts:Nt
        simulation_loop!(node_queue=node_queue, operator=operator, options=options, containers=containers, mass_flows=mass_flows, fluid_T=fluid_T, i=i)
    end
end

function first_zero_column(M)
    for j in 1:size(M, 2)
        if @views all(M[:, j] .== 0)
            return j
        end
    end
    return 1
end

"""
    simulate_steps!(;n, initial_step = nothing, options::SimulationOptions, operator, containers::SimulationContainers)

Similar to [`simulate!`](@ref), but only run `n` steps, starting at the step `initial_step`, of the simulation defined by `options`. 
"""
function simulate_steps!(;n, initial_step = nothing, operator, options::SimulationOptions, containers::SimulationContainers)
    @unpack X = containers
    @unpack medium = options
    @unpack Nb = options
    
    mass_flows = zeros(Nb+2)
    if isnothing(initial_step) || initial_step == 0
        initial_step = first_zero_column(X)
    end
    fluid_T = initial_step == 1 ? get_T0(medium) .* ones(2Nb) : X[1:2Nb, initial_step-1]

    node_queue = Queue{Tuple{Int, Float64}}()
    for i = initial_step:initial_step+n-1
        simulation_loop!(node_queue=node_queue, operator=operator, options=options, containers=containers, mass_flows=mass_flows, fluid_T=fluid_T, i=i)
    end
end

function simulation_loop!(;node_queue, operator, options::SimulationOptions, containers::SimulationContainers, mass_flows, fluid_T, i)
    @unpack configurations, method, constraint, borefield, medium, fluid, boundary_condition, approximation = options
    @unpack Nb, Ns, Nt, Ts = options
    @unpack M, b, X, mf = containers 

    operation = @views operate(operator, i, options, X[:, 1:i-1])
    operation = unwrap(operation)
    @unpack network = operation
    compute_mass_flows!(node_queue, mass_flows, network, operation)

    @views @. mf[:, i] = mass_flows[1:Nb]

    Nbr = n_branches(network)

    internal_model_eqs = 1:Nb
    topology_eqs = Nb+1:2Nb-Nbr
    constraints_eqs = 2Nb-Nbr+1:2Nb
    method_eqs = 2Nb+1:3Nb
    balance_eqs = 3Nb+1:4Nb

    # Update M
    @views internal_model_coeffs!(M[internal_model_eqs, :], borefield, medium, mass_flows, fluid_T, fluid)
    @views topology_coeffs!(M[topology_eqs, :], network, mass_flows)
    @views constraints_coeffs!(M[constraints_eqs, :], constraint, borefield, network, mass_flows)
    if i == Ts
        @views method_coeffs!(M[method_eqs, :], method, options)
    end
    @views heat_balance_coeffs!(M[balance_eqs, :], borefield, mass_flows, fluid)

    # Update b
    @views internal_model_b!(b[internal_model_eqs], borefield)
    @views topology_b!(b[topology_eqs], operation)
    @views constraints_b!(b[constraints_eqs], constraint, network, mass_flows, i)
    @views method_b!(b[method_eqs], method, borefield, medium, i)
    @views heat_balance_b!(b[balance_eqs], borefield)  

    # Solve system of equations
    @views solve_step!(X[:, i], M, b)

    # Update auxiliaries
    update_auxiliaries!(method, X, borefield, i)

    @views @. fluid_T = X[1:2Nb, i]
end
