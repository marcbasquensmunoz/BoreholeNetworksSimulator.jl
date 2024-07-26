"""
    HeatLoadConstraint(Q_tot::Matrix{T}){T <: Number} <: Constraint

Constrain the heat extracted per branch.

The heat constraint `Q_tot` must be a `Matrix`, whose column `i` are the loads per branch at the 
time step `i`. The amount of rows of `Q_tot` must equal to the amount of branches specified 
in [`BoreholeNetwork`](@ref).
"""
struct HeatLoadConstraint{T <: Number} <: Constraint
    Q_tot::Matrix{T}
end

function constant_HeatLoadConstraint(Q_tot::Vector{T}, Nt) where {T <: Number}
    Nbr = length(Q_tot)
    Q = zeros(Nbr, Nt)
    for i in 1:Nt
        Q[:, i] .= Q_tot
    end
    HeatLoadConstraint(Q)
end

function constraints_coeffs!(M, ::HeatLoadConstraint, operation)
    M .= 0
    Nb = sum([length(branch) for branch in operation.network.branches])

    for (i, branch) in enumerate(operation.network.branches)
        M[i, 3Nb .+ branch] .= 1.
    end
end

function constraints_b!(b, constraint::HeatLoadConstraint, operation, step)
    b .= constraint.Q_tot[:, step]
end

function branch_of_borehole(network, borehole)
    for (i, branch) in enumerate(network)
        if borehole in branch
            return i
        end
    end
    return 0
end