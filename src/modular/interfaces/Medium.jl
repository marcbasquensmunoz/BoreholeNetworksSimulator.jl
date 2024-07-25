
"""
    abstract type Medium

Interface for mediums.

Required functions:
- `get_λ(::Medium)`: Return the thermal conductivity of the medium.
- `get_α(::Medium)`: Return the thermal diffusivity of the medium.
- `get_T0(::Medium)`: Return the initial temperature of the medium.
- `compute_response!(g, ::Medium, borefield::Borefield, boundary_condition::BoundaryCondition, t)`: 
    Compute inplace in `g` the thermal responses between boreholes in `borefield`, 
    imposing the boundary condition `boundary_condition`, for all times in `t`.
"""
abstract type Medium end

@required Medium begin
    get_λ(::Medium)
    get_α(::Medium)
    get_T0(::Medium)
    compute_response!(g, ::Medium, borefield, boundary_condition, t) 
end
