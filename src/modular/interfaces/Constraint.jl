
"""
    abstract type Constraint

Common interface for constraints
"""
abstract type Constraint end

# Compute the coefficient matrix of the constraint equations
function branches_constraints_coeffs!(M, ::Constraint, operation) end
# Compute the independent vector of the constraint equations
function branches_constraints_b!(b, ::Constraint, operation, step) end