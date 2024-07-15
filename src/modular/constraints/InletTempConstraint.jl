
struct InletTempConstraint{T} <: Constraint 
    T_in::Vector{T}
end

function constraints_coeffs!(M, ::InletTempConstraint, operation)
    M .= 0.
    for (i, borehole) in enumerate(first_boreholes(operation.network))
        M[i, 2*borehole - 1] = 1.
    end
end

function constraints_b!(b, constraint::InletTempConstraint, operation, step)
    b .= constraint.T_in
end
