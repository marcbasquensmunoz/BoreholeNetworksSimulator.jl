
"""
    abstract type Operator

Interface for operation strategies.

Required functions:
- `operate(::Operator, i, options, Tfin, Tfout, Tb, q)`: Return a `BoreholeOperation`.
"""
abstract type Operator end

@required Operator begin
    operate(::Operator, i, options, Tfin, Tfout, Tb, q)
end
