using GeometryTypes: Point2
using Parameters
using CSV, DataFrames

@with_kw struct BoreholeOperation
    network
    mass_flows
end

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


function heat_balance_coeffs!(M, borefield::Borefield, operation::BoreholeOperation, cpf)
    Nb = borehole_amount(borefield)
    Ns = segment_amount(borefield)

    for i in 1:Nb
        M[i, i*2-1:i*2] = [cpf -cpf] .* operation.mass_flows[i]
    end

    map = segment_map(borefield)
    for i in 1:Nb
        for j in 1:Ns
            if map[j] == i
                M[i, 3Nb+j] = -get_h(borefield, i)
            end
        end
    end   
end

function heat_balance_b!(b, borefield, current_Q)
    Nb = borehole_amount(borefield)
    for i in 1:Nb
        b[i] = current_Q[i] * get_H(borefield, i)
    end  
end

function solve_step!(X, A, b, step, Nb, current_Q)
    x = A\b
    X[step,:] = x
    current_Q .+= x[3Nb+1:end]
end


function load_borefield_from_file(file)
    df = CSV.File(file; decimal=',', delim = ';') |> DataFrame
    borehole_positions =  [Point2(x, y) for (x,y) in zip(df.X,df.Y) ]
    EqualBoreholesBorefield(borehole_prototype=SingleUPipeBorehole(H=50., D=4.), positions=borehole_positions, medium=GroundWaterMedium())
end

function operator(i, Tin, Tout, Tb, Δq, Q)
    BoreholeOperation(networks[i%12 in 1:6 ? 2 : 1], ones(48))
end

networks = 
[
    [   
        [22,30,37,29,36,35], 
        [34,40,41,42,48,43],  
        [33,32,39,45,46,47],                        
        [26,25,24,31,38,44],  
        [23,16,17,11,12,18],  
        [10,5,6,7,13,19],     
        [1,2,3,8,14,20],      
        [4,9,15,21,28,27]                                         
    ]
    ,
    [
        [35,36,29,37,30,22],
        [43,48,42,41,40,34],
        [47,46,45,39,32,33],
        [44,38,31,24,25,26],
        [18,12,11,17,16,23],
        [19,13,7,6,5,10],
        [20,14,8,3,2,1],
        [27,28,21,15,9,4]
    ]
]
tstep = 8760*3600/12.
tmax  = 8760*3600*10.
Nt = div(tmax, tstep)

constraint = InletTempConstraint([i%12 in 1:6 ? 90. : 55. for i = 1:Nt])
model = ConvolutionGroundModel(T0 = 10.)
cdir = @__DIR__
borehole_positions_file = "$cdir/../../examples/example1/data/Braedstrup_borehole_coordinates.txt"
borefield = load_borefield_from_file(borehole_positions_file)

simulate(operator=operator, borefield=borefield, constraint=constraint, model=model, tstep=tstep, tmax=tmax)