
"""
    abstract type Fluid end

Interface for fluids.

Required functions:
- `cpf(::Fluid)`: Return the scpecific heat capacity of the fluid.
- `thermophysical_properties(::Fluid, T)`: Return the `μ`, `ρ`, `cp`, and `k` of the fluid at temperature `T`.
"""
abstract type Fluid end

@required Fluid begin
    cpf(::Fluid)
    thermophysical_properties(::Fluid, T)
end
