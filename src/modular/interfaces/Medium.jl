
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



"""
    MediumMock <: Medium

Mock for testing purposes.
"""
@with_kw struct MediumMock <: Medium 
    λ = 0.
    α = 0.
    T0 = 0.
    g = 0.
    q_coef = 0.
end
get_λ(m::MediumMock) = m.λ
get_α(m::MediumMock) = m.α
get_T0(m::MediumMock) = m.T0
compute_response!(g, m::MediumMock, borefield, boundary_condition, t) = g .= m.g
constant_integral(::MediumMock, method, setup, λ, i) = 0.