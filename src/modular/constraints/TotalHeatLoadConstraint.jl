"""
    TotalHeatLoadConstraint(Q_tot::Vector{T}){T <: Number} <: Constraint

Constrain the total heat extracted in the borefiled.

The heat constraint `Q_tot` must be a `Vector`, whose elements are the total load at the time step `i`.
"""
struct TotalHeatLoadConstraint{T <: Number} <: Constraint
    Q_tot::Vector{T}
end


function constraints_coeffs!(M, ::TotalHeatLoadConstraint, borefield::Borefield, network, mass_flows)
    M .= zero(eltype(M))
    Nb = n_boreholes(network)
    j = 1

    first_bhs = first_bhs_in_branch(network)
    first_bh = findfirst(i-> mass_flows[i] != 0., first_bhs)

    for bh in first_bhs
        if mass_flows[bh] == 0.
            M[j, 2*bh] = -1.
        else 
            if bh == first_bh
                continue
            end
            M[j, 2*first_bh-1] = -1.
        end
        M[j, 2*bh-1] = 1.
        j += 1
    end

    if sum(mass_flows) != 0
        for i in 1:Nb
            M[end, 3Nb+i] = get_H(borefield, i)
        end
    end
end

function constraints_b!(b, constraint::TotalHeatLoadConstraint, network, mass_flows, step)
    b .= zero(eltype(b))
    if sum(mass_flows) != 0.
        b[end] = constraint.Q_tot[step]
    end
end
