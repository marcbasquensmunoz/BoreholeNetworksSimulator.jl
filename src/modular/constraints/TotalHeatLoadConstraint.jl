"""
    TotalHeatLoadConstraint(Q_tot::Vector{T}){T <: Number} <: Constraint

Constrain the total heat extracted in the borefiled.

The heat constraint `Q_tot` must be a `Vector`, whose elements are the total load at the time step `i`.
"""
struct TotalHeatLoadConstraint{T <: Number} <: Constraint
    Q_tot::Vector{T}
end


function constraints_coeffs!(M, ::TotalHeatLoadConstraint, operation, borefield)
    M .= zero(eltype(M))

    Nb = n_boreholes(operation.network)

    first_bh = operation.network.branches[1][1]
    for (i, branch) in enumerate(operation.network.branches[2:end])
        M[i+1, 2*first_bh-1] = -1.
        bh_in = 2*branch[1]-1
        M[i+1, bh_in] = 1.
    end

    first_bh = 4
    if sum(abs.(operation.mass_flows)) == 0
        M[1, 2*first_bh-1] = 1.
        M[1, 2Nb+first_bh] = -1.
        return
    end

    M[1, 3Nb+1:end] .= [get_H(borefield, i) for i in 1:Nb]
end

function constraints_b!(b, constraint::TotalHeatLoadConstraint, operation, step)
    b .= zero(eltype(b))
    b[1] = constraint.Q_tot[step]
end
