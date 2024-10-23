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

    first_bhs = first_bhs_in_branch(network)
    first_bh = first_bhs[1]
    for i in 2:n_branches(network)
        M[i, 2*first_bh-1] = -1.
        bh_in = 2*first_bhs[i]-1
        M[i, bh_in] = 1.
    end

    if sum(mass_flows) == 0
        M[1, 2*first_bh-1] = 1.
        M[1, 2Nb+first_bh] = -1.
        return
    end

    for i in 1:Nb
        M[1, 3Nb+i] = get_H(borefield, i)
    end
end

function constraints_b!(b, constraint::TotalHeatLoadConstraint, operation, step)
    b .= zero(eltype(b))
    b[1] = constraint.Q_tot[step]
end
