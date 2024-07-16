struct HeatLoadConstraint{T} <: Constraint
    Q_tot::Vector{T}
end

function constraints_coeffs!(M, ::HeatLoadConstraint, operation)
    M .= 0
    Nb = sum([length(branch) for branch in operation.network.branches])

    for (i, branch) in enumerate(operation.network.branches)
        M[i, 3Nb .+ branch] .= 1.
    end

    #=
    first_branch = operation.network[1]
    for i = 1:2*Nb
        branch = branch_of_borehole(operation.network, div(i-1, 2) + 1)
        M[first_branch[1], i] = (-1)^i * operation.mass_flows[branch] * operation.cpf
    end

    for i in 2:Nbr   
        branch = operation.network[i]
        M[branch[1], 2*first_branch[1] - 1] = 1.
        M[branch[1], 2*branch[1] - 1] = -1.
    end
    for branch in operation.network     
        for i in 2:length(branch)
            M[branch[i], 2*branch[i] - 1] = 1.
            M[branch[i], 2*branch[i-1]] = -1.
        end
    end
    M[1, 1] = 0.
    M[1, 2] = 0.
    M[1, 4] = 1.
    =#
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