
"""
    EqualBoreholesBorefield{T <: Borehole, S <: Number} <: Borefield
    EqualBoreholesBorefield(borehole_prototype::T, positions::Vector{Point2{S}}))

Model a borefield with boreholes all identical to the prototype `borehole_prototype`, placed at `positions`.
Note that the length of `positions` determines the amount of boreholes in the field.
"""
@with_kw struct EqualBoreholesBorefield{T <: Borehole, S <: Number} <: Borefield
    borehole_prototype::T
    positions::Vector{Tuple{S, S}}
    Nb::Int = length(positions)
    initial_ΔT::Vector{S} = zeros(S, Nb)
end

n_boreholes(bf::EqualBoreholesBorefield) = bf.Nb
n_segments(bf::EqualBoreholesBorefield) = bf.Nb
get_H(bf::EqualBoreholesBorefield, i) = get_H(bf.borehole_prototype)
get_h(bf::EqualBoreholesBorefield, i) = get_h(bf.borehole_prototype)
get_rb(bf::EqualBoreholesBorefield, i) = get_rb(bf.borehole_prototype)
where_is_segment(bf::EqualBoreholesBorefield, i) = div((i-1), get_n_segments(bf.borehole_prototype)) + 1  
T_past_influence(bf::EqualBoreholesBorefield) = bf.initial_ΔT

function segment_coordinates(bf::EqualBoreholesBorefield, segment)
    D = get_D(bf.borehole_prototype)
    h = get_h(bf.borehole_prototype)

    position = bf.positions[where_is_segment(bf, segment)]
    i = mod((segment-1), get_n_segments(bf.borehole_prototype)) + 1
    z_ref = D + (i - 1) * h
    
    (position[1], position[2], z_ref, h)
end

function internal_model_coeffs!(M, borefield::EqualBoreholesBorefield, medium, mass_flows, T_fluid, fluid)
    Nb = n_boreholes(borefield)
    λ = get_λ(medium)

    for i in 1:Nb
        mass_flow = mass_flows[i]

        Tref = (T_fluid[2i - 1] + T_fluid[2i]) / 2

        k_in, k_out, k_b = uniform_Tb_coeffs(borefield.borehole_prototype, λ, mass_flow, Tref, fluid)

        M[i, 2i - 1]  = k_in
        M[i, 2i]      = k_out
        M[i, 2Nb + i] = k_b    
    end
end

function internal_model_b!(b, ::EqualBoreholesBorefield) end
