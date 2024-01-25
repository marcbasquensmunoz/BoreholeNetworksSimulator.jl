struct InletTempConstraint <: Constraint 
    T_in
end

function branches_constraints_coeffs!(M, ::InletTempConstraint, operation)
    M .= 0
    Nb = sum([length(branch) for branch in operation.network])

    for i = 1:Nb
        M[i, 2i-1] = 1
    end
    for branch in operation.network     
        for i in 2:length(branch)
            M[branch[i], (branch[i-1])*2] = -1 
        end
    end
end

function branches_constraints_b!(b, constraint::InletTempConstraint, operation, step)
    b .= 0
    for branch in operation.network              
        b[branch[1]] = constraint.T_in[step]
    end
end
