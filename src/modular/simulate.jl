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

`operator` should be a function that returns a `BoreholeOperation` and with signature `operator(i, Tin, Tout, Tb, q)`:
- `i::Int` is the time step
- `Tin` is a vector containing the inlet temperature of each borehole
- `Tout` is a vector containing the outlet temperature of each borehole
- `Tb` is a vector containing the borehole wall temperature of each borehole
- `q` is a vector containing the heat exchanged by each borehole
"""
function simulate!(;operator, options::SimulationOptions, containers::SimulationContainers)
    @unpack method, constraint, borefield, medium, fluid, boundary_condition = options
    @unpack Nb, Ns, Nt, Ts = options
    @unpack M, b, X = containers 

    compatibility = check_compatibility(medium, constraint, method)
    if compatibility isa NotCompatible
        println(compatibility.message)
        return
    end
    
    last_operation = BoreholeOperation(nothing)

    # Simulation loop
    for i = Ts:Nt
        operation = @views operator(i, X[1:2:2Nb, 1:i], X[2:2:2Nb, 1:i], X[2Nb+1:3Nb, 1:i], X[3Nb+1:end, 1:i])
        if ispy(operation)
            operation = PythonCall.pyconvert(BoreholeOperation, operation)
        end

        Nbr = n_branches(operation.network)

        internal_model_eqs = 1:Nb
        topology_eqs = Nb+1:2Nb-Nbr
        constraints_eqs = 2Nb-Nbr+1:2Nb
        method_eqs = 2Nb+1:3Nb
        balance_eqs = 3Nb+1:4Nb

        # Update M
        @views internal_model_coeffs!(M[internal_model_eqs, :], borefield, medium, operation, i == 1 ? get_T0(medium) .* ones(2Nb) :  X[1:2Nb, i-1], fluid)
        if last_operation.network != operation.network
            @views topology_coeffs!(M[topology_eqs, :], operation)
        end
        @views constraints_coeffs!(M[constraints_eqs, :], constraint, operation)
        if i == Ts
            @views method_coeffs!(M[method_eqs, :], method, borefield, medium, boundary_condition)
        end
        if last_operation.mass_flows != operation.mass_flows
            @views heat_balance_coeffs!(M[balance_eqs, :], borefield, operation, fluid)
        end

        # Update b
        @views internal_model_b!(b[internal_model_eqs], borefield)
        @views constraints_b!(b[constraints_eqs], constraint, operation, i)
        @views method_b!(b[method_eqs], method, borefield, medium, i)
        @views heat_balance_b!(b[balance_eqs], borefield, X[balance_eqs])  

        # Solve system of equations
        @views solve_step!(X[:, i], M, b)

        # Update auxiliaries
        update_auxiliaries!(method, X, borefield, i)

        last_operation = operation
    end
end
