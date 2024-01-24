using Parameters

@with_kw struct EqualBoreholesBorefield <: Borefield
    borehole_prototype::Borehole
    positions
    medium::Medium
end

borehole_amount(bf::EqualBoreholesBorefield) = length(bf.positions)
segment_amount(bf::EqualBoreholesBorefield) = length(bf.positions)
get_H(bf::EqualBoreholesBorefield, i) = get_H(bf.borehole_prototype)
get_h(bf::EqualBoreholesBorefield, i) = get_h(bf.borehole_prototype)
get_rb(bf::EqualBoreholesBorefield, i) = get_rb(bf.borehole_prototype)
segment_map(bf::EqualBoreholesBorefield) = 1:length(bf.positions)           # Return an array indicating to which borehole each segments belongs

function segment_coordinates(bf::EqualBoreholesBorefield)
    D = get_D(bf.borehole_prototype)
    H = get_H(bf.borehole_prototype)
    h = get_h(bf.borehole_prototype)

    z_ref = collect(D:h:D+H-h)          # line source reference point
    z_eval = collect(D+h/2:h:D+H-h/2)   # evaluation points (evaluate at the mid point of the segment)
      
    coord_source = [(x[1],x[2],p) for x in bf.positions for p in z_ref]   # position of sources 
    coord_eval = [(x[1],x[2],p) for x in bf.positions for p in z_eval]    # position of evaluation points

    return (coord_source, coord_eval)
end

function internal_model_coeffs!(M, borefield::EqualBoreholesBorefield, operation, cpf)
    Nb = borehole_amount(borefield)
    for i = 1:Nb
        mass_flow = operation.mass_flows[i]
        R = resistance_network(borefield.borehole_prototype, λ(borefield.medium), mass_flow)
        A = coefficient_matrix(R, cpf, mass_flow)
        k_in, k_out, k_b = uniformTb_koeff(A, get_H(borefield.borehole_prototype)) 

        M[i, i*2 - 1]  = k_in[1]
        M[i, i*2]      = k_out[1]
        M[i, Nb*2 + i] = k_b[1]    
    end
end

function internal_model_b!(b, borefield::EqualBoreholesBorefield) end
