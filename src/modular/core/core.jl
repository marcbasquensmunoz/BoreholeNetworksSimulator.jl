using Parameters

@with_kw struct BoreholeOperation
    network
    mass_flows
    cpf = 4182.        # specific heat capacity
end
BoreholeOperation(::Nothing) =  BoreholeOperation(nothing, nothing, nothing)

function heat_balance_coeffs!(M, borefield::Borefield, operation::BoreholeOperation)
    Nb = borehole_amount(borefield)
    Ns = segment_amount(borefield)

    for i in 1:Nb
        M[i, i*2-1:i*2] = operation.cpf .* operation.mass_flows[i] .* [1 -1]
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
