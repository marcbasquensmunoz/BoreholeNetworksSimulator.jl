"""
    HeatLoadConstraint(Q_tot::Vector) <: Constraint

Constrain to `Q_tot` the heat extracted per branch.
Note that `Q_tot` must have length equal to the amount of branches specified in [`BoreholeNetwork`](@ref).
"""
struct HeatLoadConstraint <: Constraint
    Q_tot::Vector
end

function constraints_coeffs!(M, ::HeatLoadConstraint, operation)
    M .= 0
    Nb = sum([length(branch) for branch in operation.network.branches])

    for (i, branch) in enumerate(operation.network.branches)
        M[i, 3Nb .+ branch] .= 1.
    end
end

function constraints_b!(b, constraint::HeatLoadConstraint, operation, step)
    b .= constraint.Q_tot
end

function branch_of_borehole(network, borehole)
    for (i, branch) in enumerate(network)
        if borehole in branch
            return i
        end
    end
    return 0
end