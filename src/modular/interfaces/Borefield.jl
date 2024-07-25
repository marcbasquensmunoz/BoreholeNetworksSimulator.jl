"""
    abstract type Borefield

Interface for borefields.

Required functions
- `n_boreholes(::Borefield)`: Return the amount of boreholes present in the borefield.
- `get_H(::Borefield, i)`: Return the length of borehole `i`.
- `get_rb(::Borefield, i)`: Return the radius of borehole `i`.
- `segment_coordinates(::Borefield)`: Return a vector with the coordinates of each segment.  
- `internal_model_coeffs!(M, ::Borefield, medium, operation, fluid)`: Compute inplace in `M`
    the coefficients corresponding to the internal model equations, given the `medium`, `fluid` and
    `operation` in use in the system.
    Note that `M` is only a slice of `Nb` (number of boreholes) rows, provided as a `view`.
- `internal_model_b!(b, ::Borefield)`: Compute inplace in `b`
    the independent terms corresponding to the internal model equations.
    Note that `b` is only a vector of length `Nb` (number of boreholes) rows, provided as a `view`.
"""
abstract type Borefield end

@required Borefield begin
    n_boreholes(::Borefield)
    get_H(::Borefield, i)
    get_rb(::Borefield, i)
    segment_coordinates(::Borefield, segment)
    internal_model_coeffs!(M, ::Borefield, medium, operation, T_fluid, fluid)
    internal_model_b!(b, ::Borefield)
end
