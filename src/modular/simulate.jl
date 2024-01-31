using JLD2

function simulate(;operator, parameters::SimulationParameters, containers::SimulationContainers, borefield::Borefield, constraint::Constraint, method::Method)

    @unpack Nb, Ns, Nt, Ts = parameters
    @unpack M, b, X, current_Q = containers 
    
    last_operation = BoreholeOperation(nothing)

    # Simulation loop
    for i = Ts:Nt
        operation = @views operator(i, X[:, 1:2:2Nb], X[:, 2:2:2Nb], X[:, 2Nb+1:3Nb], X[3Nb+1:end], current_Q)

        # Update M
        if last_operation.mass_flows != operation.mass_flows
            @views internal_model_coeffs!(M[1:Nb, :], borefield, operation)
        end
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
        @views method_b!(b[2Nb+1:2Nb+Ns], method, borefield, i)
        @views heat_balance_b!(b[2Nb+Ns+1:3Nb+Ns], borefield, current_Q)  

        # Solve system of equations
        solve_step!(X, M, b, i, Nb, current_Q)

        # Update auxiliaries
        update_auxiliaries!(method, X, borefield, i)

        last_operation = operation
    end
end
