
@with_kw struct EqualBoreholesBorefield{T <: Borehole, R <: Medium, S <: Real} <: Borefield
    borehole_prototype::T
    positions::Vector{Point2{Float64}}
    medium::R
    T0::S
end

borehole_amount(bf::EqualBoreholesBorefield) = length(bf.positions)
segment_amount(bf::EqualBoreholesBorefield) = length(bf.positions)
get_H(bf::EqualBoreholesBorefield, i) = get_H(bf.borehole_prototype)
get_h(bf::EqualBoreholesBorefield, i) = get_h(bf.borehole_prototype)
get_rb(bf::EqualBoreholesBorefield, i) = get_rb(bf.borehole_prototype)
get_T0(bf::EqualBoreholesBorefield) = bf.T0
where_is_segment(bf::EqualBoreholesBorefield, i) = div((i-1), get_n_segments(bf.borehole_prototype)) + 1  

function segment_coordinates(bf::EqualBoreholesBorefield, segment)
    D = get_D(bf.borehole_prototype)
    h = get_h(bf.borehole_prototype)

    position = bf.positions[where_is_segment(bf, segment)]
    i = mod((segment-1), get_n_segments(bf.borehole_prototype)) + 1
    z_ref = D + (i - 1) * h
    z_eval = z_ref + h/2

    # (x, y, z_ref, z_eval)
    return (position[1], position[2], z_ref, z_eval)
end

function internal_model_coeffs!(M, borefield::EqualBoreholesBorefield, operation, T_fluid)
    Nb = borehole_amount(borefield)

    for (i, branch) in enumerate(operation.network)
        mass_flow = operation.mass_flows[i]

        for j in branch
            Tref = (T_fluid[2*j - 1] + T_fluid[2*j]) / 2

            k_in, k_out, k_b = uniform_Tb_coeffs(borefield.borehole_prototype, get_λ(borefield.medium), mass_flow, Tref, operation.cpf)

            M[j, j*2 - 1]  = k_in
            M[j, j*2]      = k_out
            M[j, Nb*2 + j] = k_b    
        end
    end
end

function internal_model_b!(b, borefield::EqualBoreholesBorefield) end
