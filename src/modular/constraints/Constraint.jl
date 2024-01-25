"""
branches_constraints_coeffs!(M, constraint, operation)                      Compute the coefficient matrix of the constraint equations
function branches_constraints_b!(b, constraint, operation, step)            Compute the independent vector of the constraint equations
"""
abstract type Constraint end