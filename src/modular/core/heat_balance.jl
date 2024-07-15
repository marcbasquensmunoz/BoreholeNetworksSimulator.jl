
function heat_balance_coeffs!(M, borefield::Borefield, operation::BoreholeOperation)
    Nb = borehole_amount(borefield)
    Ns = segment_amount(borefield)

    for i in 1:Nb
        heat = operation.cpf .* operation.mass_flows[branch_of_borehole(operation, i)] 
        M[i, i*2-1] = heat
        M[i, i*2]  = -heat
    end

    for i in 1:Nb
        for j in 1:Ns
            if where_is_segment(borefield, j) == i
                M[i, 3Nb+j] = -get_h(borefield, i)
            end
        end
    end   
end

function heat_balance_b!(b, borefield, Q)
    Nb = borehole_amount(borefield)
    for i in 1:Nb
        b[i] = Q[i] * get_H(borefield, i)
    end  
end
