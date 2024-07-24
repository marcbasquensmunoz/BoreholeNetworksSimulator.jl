
"""
    EqualBoreholesBorefield{T <: Borehole, R <: Medium, S <: Real} <: Borefield
    EqualBoreholesBorefield(borehole_prototype::T, positions::Vector{Point2{S}}), medium::R)

Model a borefield with boreholes all identical to the prototype `borehole_prototype`, placed at `positions`.
Note that the length of `positions` determines the amount of boreholes in the field.
`medium` contains the properties of the ground.
"""
@with_kw struct EqualBoreholesBorefield{T <: Borehole, S <: Real} <: Borefield
    borehole_prototype::T
    positions::Vector{Tuple{S, S}}
end

n_boreholes(bf::EqualBoreholesBorefield) = length(bf.positions)
n_segments(bf::EqualBoreholesBorefield) = length(bf.positions)
get_H(bf::EqualBoreholesBorefield, i) = get_H(bf.borehole_prototype)
get_h(bf::EqualBoreholesBorefield, i) = get_h(bf.borehole_prototype)
get_rb(bf::EqualBoreholesBorefield, i) = get_rb(bf.borehole_prototype)
where_is_segment(bf::EqualBoreholesBorefield, i) = div((i-1), get_n_segments(bf.borehole_prototype)) + 1  

function segment_coordinates(bf::EqualBoreholesBorefield, segment)
    D = get_D(bf.borehole_prototype)
    h = get_h(bf.borehole_prototype)

    position = bf.positions[where_is_segment(bf, segment)]
    i = mod((segment-1), get_n_segments(bf.borehole_prototype)) + 1
    z_ref = D + (i - 1) * h
    z_eval = z_ref + h/2

    # (x, y, z_ref, z_eval)
    #return (position[1], position[2], z_ref, z_eval)
    (position[1], position[2], z_ref, h)
end

function internal_model_coeffs!(M, borefield::EqualBoreholesBorefield, medium::Medium, operation, T_fluid, fluid)
    Nb = n_boreholes(borefield)

    for (i, branch) in enumerate(operation.network.branches)
        mass_flow = operation.mass_flows[i]

        for j in branch
            Tref = (T_fluid[2*j - 1] + T_fluid[2*j]) / 2

            k_in, k_out, k_b = uniform_Tb_coeffs(borefield.borehole_prototype, get_λ(medium), mass_flow, Tref, fluid)

            M[j, j*2 - 1]  = k_in
            M[j, j*2]      = k_out
            M[j, Nb*2 + j] = k_b    
        end
    end
end

function internal_model_b!(b, ::EqualBoreholesBorefield) end
