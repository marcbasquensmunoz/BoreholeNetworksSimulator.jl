
"""
    abstract type Medium

Interface for mediums.

Required functions:
- `get_λ(::Medium)`: Return the thermal conductivity of the medium.
- `get_α(::Medium)`: Return the thermal diffusivity of the medium.
- `get_T0(::Medium)`: Return the initial temperature of the medium.
"""
abstract type Medium end

@required Medium begin
    get_λ(::Medium)
    get_α(::Medium)
    get_T0(::Medium)
end
