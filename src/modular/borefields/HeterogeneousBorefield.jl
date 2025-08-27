
@with_kw struct HeterogeneousBorefield{S <: Number} <: Borefield
    boreholes::Vector{Borehole}
    positions::Vector{Tuple{S, S}}
    Nb::Int = length(positions)
end

n_boreholes(bf::HeterogeneousBorefield) = bf.Nb
get_H(bf::HeterogeneousBorefield, i) = get_H(bf.boreholes[i])
get_h(bf::HeterogeneousBorefield, i) = get_h(bf.boreholes[i])
get_rb(bf::HeterogeneousBorefield, i) = get_rb(bf.boreholes[i])

function segment_coordinates(bf::HeterogeneousBorefield, segment)
    @unpack boreholes, Nb = bf
    borehole = boreholes[segment]

    D = get_D(borehole)
    h = get_h(borehole)

    position = bf.positions[segment]
    i = mod((segment-1), get_n_segments(borehole)) + 1
    z_ref = D + (i - 1) * h
   
    (position[1], position[2], z_ref, h)
end

function internal_model_coeffs!(M, borefield::HeterogeneousBorefield, medium, mass_flows, T_fluid, fluid)
    @unpack boreholes, Nb = borefield
    λ = get_λ(medium)

    for i in 1:Nb
        mass_flow = mass_flows[i]

        Tref = (T_fluid[2i - 1] + T_fluid[2i]) / 2

        k_in, k_out, k_b = uniform_Tb_coeffs(boreholes[i], λ, mass_flow, Tref, fluid)

        M[i, 2i - 1]  = k_in
        M[i, 2i]      = k_out
        M[i, 2Nb + i] = k_b    
    end
end

function internal_model_b!(b, ::HeterogeneousBorefield) end
