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

`operator` should be a function that returns a `BoreholeOperation` and with signature `operator(i, Tin, Tout, Tb, q, configurations)`:
- `i::Int` is the time step
- `Tin` is a vector containing the inlet temperature of each borehole
- `Tout` is a vector containing the outlet temperature of each borehole
- `Tb` is a vector containing the borehole wall temperature of each borehole
- `q` is a vector containing the heat exchanged by each borehole
- `configurations`: is the list of possible hydraulic configurations of the borefield.
"""
function simulate!(;operator, options::SimulationOptions, containers::SimulationContainers)
    @unpack configurations, method, constraint, borefield, medium, fluid, boundary_condition, approximation = options
    @unpack Nb, Ns, Nt, Ts = options
    @unpack M, b, X, mf = containers 
    
    last_operation = BoreholeOperation(nothing)
    mass_flows = zeros(Nb+2)
    fluid_T = get_T0(medium) .* ones(2Nb)

    # Simulation loop
    for i = Ts:Nt
        operation = @views operate(operator, i, options, X[:, 1:i-1])
        operation = unwrap(operation)
        @unpack network = operation
        compute_mass_flows!(mass_flows, network, operation)

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
        #if last_operation.mass_flows != operation.mass_flows
            @views heat_balance_coeffs!(M[balance_eqs, :], borefield, mass_flows, fluid)
        #end

        # Update b
        @views internal_model_b!(b[internal_model_eqs], borefield)
        @views topology_b!(b[topology_eqs], operation)
        @views constraints_b!(b[constraints_eqs], constraint, operation, i)
        @views method_b!(b[method_eqs], method, borefield, medium, i)
        @views heat_balance_b!(b[balance_eqs], borefield)  

        # Solve system of equations
        @views solve_step!(X[:, i], M, b)

        # Update auxiliaries
        update_auxiliaries!(method, X, borefield, i)

        @views @. fluid_T = X[1:2Nb, i]
        last_operation = operation
    end
end
