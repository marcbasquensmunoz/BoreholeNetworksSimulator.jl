using JLD2

function simulate(;parameters::SimulationParameters, containers::SimulationContainers, operator, borefield::Borefield, constraint::Constraint, model::GroundModel, symtitle="simulation")
    
    @unpack Nb, Ns, Nt, t, Ts = parameters
    @unpack M, b, X, current_Q = containers 
    
    precompute_auxiliaries!(model, borefield, t)

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
        @views ground_model_coeffs!(M[2Nb+1:2Nb+Ns, :], model, borefield)
        if last_operation.mass_flows != operation.mass_flows
            @views heat_balance_coeffs!(M[2Nb+Ns+1:3Nb+Ns, :], borefield, operation)
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

    #=
    cdir = @__DIR__
    results_directory = "$cdir/results"
    simulation_results_directory = "$results_directory/$symtitle"
    !isdir(results_directory) && mkdir(results_directory)
    !isdir(simulation_results_directory) && mkdir(simulation_results_directory)
    save("$(simulation_results_directory)/cache_$(tmax).jld2" , 
        Dict( 
            "X" => X,
            "b" => b,
            "current_Q" => current_Q
        )
    )
    =#
end
