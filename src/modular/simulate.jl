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

function simulate(;operator, parameters::SimulationParameters, containers::SimulationContainers, borefield::Borefield, constraint::Constraint, method::Method)

    @unpack Nb, Ns, Nt, Ts = parameters
    @unpack M, b, X = containers 

    compatibility = check_compatibility(borefield, constraint, method)
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

        # Update M
        @views internal_model_coeffs!(M[1:Nb, :], borefield, operation, i == 1 ? get_T0(borefield) .* ones(2Nb) :  X[1:2Nb, i-1])
        if last_operation.network != operation.network
            @views branches_constraints_coeffs!(M[Nb+1:2Nb, :], constraint, operation)
        end
        if i == Ts
            @views method_coeffs!(M[2Nb+1:2Nb+Ns, :], method, borefield)
        end
        if last_operation.mass_flows != operation.mass_flows
            @views heat_balance_coeffs!(M[2Nb+Ns+1:3Nb+Ns, :], borefield, operation)
        end

        # Update b
        @views internal_model_b!(b[1:Nb], borefield)
        @views branches_constraints_b!(b[Nb+1:2Nb], constraint, operation, i)
        @views method_b!(b[2Nb+1:2Nb+Ns], method, borefield, i, X[3Nb+1:3Nb+Ns])
        @views heat_balance_b!(b[2Nb+Ns+1:3Nb+Ns], borefield, X[3Nb+1:3Nb+Ns])  

        # Solve system of equations
        solve_step!(X, M, b, i, Nb)

        # Update auxiliaries
        update_auxiliaries!(method, X, borefield, i)

        last_operation = operation
    end
end
