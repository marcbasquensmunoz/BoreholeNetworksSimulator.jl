"""
    InletTempConstraint(T_in) <: Constraint

Constrain to `T_in` the inlet temperature of the first borehole in each branch.
Note that `T_in` must have length equal to the amount of branches specified in [`BoreholeNetwork`](@ref).
"""
struct InletTempConstraint <: Constraint 
    T_in
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
