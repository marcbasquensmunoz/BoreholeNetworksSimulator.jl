
"""
    abstract type Constraint

Interface for constraints.

Required functions:
- `constraints_coeffs!(M, ::Constraint, operation::BoreholeOperation)`: Compute inplace in `M`
    the coefficients corresponding to the constraints equations, given the current `operation.network`. 
    Note that `M` is only a slice of `Nbr` (number of branches) rows, provided as a `view`.
- `constraints_b!(b, ::Constraint, ::BoreholeOperation, step)`: Compute inplace in `b`
    the independent term corresponding to the constraints equations, given the current `operation.network`,
    at the time step `step`. 
    Note that `b` is only a vector of length `Nbr` (number of branches), provided as a `view`.
"""
abstract type Constraint end

@required Constraint begin
    constraints_coeffs!(M, ::Constraint, operation, borefield)
    constraints_b!(b, ::Constraint, operation, step)
end
