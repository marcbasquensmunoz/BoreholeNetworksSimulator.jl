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

"""
    constant_HeatLoadConstraint(Q_tot::Vector{T}, Nt) where {T <: Number}

Convenience initializer for `HeatLoadConstraint`. It creates a constant heat load constraint through all the `Nt` time steps,
where `Q_tot` are the heat load for each branch.
"""
function constant_HeatLoadConstraint(Q_tot::Vector{T}, Nt) where {T <: Number}
    Nbr = length(Q_tot)
    Q = zeros(Nbr, Nt)
    for i in 1:Nt
        Q[:, i] .= Q_tot
    end
    HeatLoadConstraint(Q)
end


"""
    uniform_HeatLoadConstraint(Q_tot::Vector{T}, Nbr) where {T <: Number}

Convenience initializer for `HeatLoadConstraint`. It creates a uniform heat load constraint along all branches,
where `T_in` are the heat load for each time step.
"""
function uniform_HeatLoadConstraint(Q_tot::Vector{T}, Nbr) where {T <: Number}
    Nt = length(Q_tot)
    Q = zeros(Nbr, Nt)
    for i in 1:Nbr
        Q[i, :] .= Q_tot
    end
    HeatLoadConstraint(Q)
end

function constraints_coeffs!(M, ::HeatLoadConstraint, borefield::Borefield, network, mass_flows)
    M .= zero(eltype(M))

    Nb = n_boreholes(network)
    for (i, first_bh) in enumerate(first_bhs_in_branch(network))
        for bh in neighborhood(network.graph, first_bh, Nb)
            if bh == source(network) || bh == sink(network) 
                continue
            end
            M[i, 3Nb + bh] = one(eltype(M)) * get_H(borefield, bh)
        end
    end
end

function constraints_b!(b, constraint::HeatLoadConstraint, operation, step)
    b .= @view constraint.Q_tot[:, step]
end

function branch_of_borehole(network, borehole)
    for (i, branch) in enumerate(network)
        if borehole in branch
            return i
        end
    end
    return 0
end