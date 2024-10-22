
"""
    abstract type Operator

Interface for operation strategies.

Required functions:
- `operate(::Operator, step, options, X)` 

The implementation of `operate` must return an instance of `BoreholeOperation`, specifying the network topology and the mass flows per branch for the current time step. 
`step` is the current time step, `options` contains the provided simulation options and `X` contains the time series up to the time step `step-1` of 
the inlet fluid temperature, the outlet fluid temperature, the wall temperature and the heat extraction rate, for each borehole in the borefield.
To get `Tfin`, `Tfout`, `Tb`, and `q` from `X`, use the functions `extract_Tfin`, `extract_Tfout`, `extract_Tb`, `extract_q`, or get them manually at the 
indices [1:2:2Nb, :], [2:2:2Nb, :], [2Nb+1:3Nb, :], [3Nb+1:end, :], respectively.
"""
abstract type Operator end

@required Operator begin    
    operate(::Operator, step, options, X)
end
