"""
    abstract type TimeSuperpositionMethod

Interface for time superposition methods.

Required functions:
- `method_coeffs!(M, ::TimeSuperpositionMethod, borefield, medium, boundary_condition)`: Compute inplace in `M`
    the coefficients corresponding to the heat transfer equations, given the `medium`, 
    and `boundary_condition` in use in the system.
    Note that `M` is only a slice of `Nbr` (number of branches) rows, provided as a `view`.
- `method_b!(b, ::TimeSuperpositionMethod, borefield, medium, step)`: Compute inplace in `b`
    the independent terms corresponding to the heat transfer equations, given the `medium`, 
    at the given time step `step`.
    Note that `b` is only a vector of length `Nbr` (number of branches) rows, provided as a `view`.
- `precompute_auxiliaries!(method::TimeSuperpositionMethod; options)`: Compute inplace in `method` 
    the auxiliary quantities used in the simulation that can be performed ahead of time.
- `update_auxiliaries!(::TimeSuperpositionMethod, X, borefield, step)`: Update inplace in `method`
    the auxiliaries after each time step `step`.
"""
abstract type TimeSuperpositionMethod end

@required TimeSuperpositionMethod begin
    method_coeffs!(M, ::TimeSuperpositionMethod, borefield, medium, boundary_condition)
    method_b!(b, ::TimeSuperpositionMethod, borefield, medium, step)
    precompute_auxiliaries!(method::TimeSuperpositionMethod, options)
    update_auxiliaries!(method::TimeSuperpositionMethod, X, borefield, step)
end
