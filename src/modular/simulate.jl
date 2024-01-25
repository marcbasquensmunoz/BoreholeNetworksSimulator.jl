
function simulate(;operator, borefield::Borefield, constraint::Constraint, model::GroundModel, tstep, tmax)

    t = tstep:tstep:tmax
    Nt = length(t)

    Nb = borehole_amount(borefield)
    Ns = segment_amount(borefield)

    cpf = 4182.        # specific heat capacity

    # Initialize system of equations
    M = zeros(3Nb+Ns, 3Nb+Ns)    
    b = zeros(3Nb + Ns)            
    X = zeros(Nt, 3Nb + Ns)          
    # Net heat injection on a given segment (this variable is needed by the solver)
    current_Q = zeros(Ns)          
    
    precompute_auxiliaries!(model, borefield, t)

    last_operation = BoreholeOperation(nothing, nothing)

    # Simulation loop
    for i = 1:Nt
        operation = @views operator(i, X[:, 1:2:2Nb], X[:, 2:2:2Nb], X[:, 2Nb+1:3Nb], X[3Nb+1:end], current_Q)

        # Add smart updating of the elements depending on what changed in the operation

        # Update M
        if last_operation.mass_flows != operation.mass_flows
            @views internal_model_coeffs!(M[1:Nb, :], borefield, operation, cpf)
        end
        if last_operation.network != operation.network
            @views branches_constraints_coeffs!(M[Nb+1:2Nb, :], constraint, operation)
        end
        @views ground_model_coeffs!(M[2Nb+1:2Nb+Ns, :], model, borefield)
        if last_operation.mass_flows != operation.mass_flows
            @views heat_balance_coeffs!(M[2Nb+Ns+1:3Nb+Ns, :], borefield, operation, cpf)
        end

        # Update b
        @views internal_model_b!(b[1:Nb], borefield)
        @views branches_constraints_b!(b[Nb+1:2Nb], constraint, operation, i)
        @views ground_model_b!(b[2Nb+1:2Nb+Ns], model, borefield, i)
        @views heat_balance_b!(b[2Nb+Ns+1:3Nb+Ns], borefield, current_Q)  

        # Solve system of equations
        solve_step!(X, M, b, i, Nb, current_Q)

        # Update auxiliaries
        update_auxiliaries!(model, X, borefield, i)

        last_operation = operation
    end

    save("$(results_directory)/cache$(symtitle)_$(tmax).jld2" , 
        Dict( 
            "X" => X,
            "b" => b,
            "current_Q" => current_Q
        )
    )
end
