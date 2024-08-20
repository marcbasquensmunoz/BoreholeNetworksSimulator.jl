"""
    GroundMedium{T <: Real} <: Medium @deftype T

Model pure conduction in the ground.

# Arguments
- `λ = 3.`: ground conductivity
- `α = 1e-6`: ground thermal diffusivity
- `C = λ/α`: ground medium capacity
- `T0 = 0.`: initial ground temperature
"""
@with_kw struct GroundMedium{T <: Real} <: Medium @deftype T
    λ = 3.               
    α = 1e-6
    T0 = 0.
    C = λ/α
end
get_λ(medium::GroundMedium) = medium.λ
get_α(medium::GroundMedium) = medium.α
get_T0(medium::GroundMedium) = medium.T0
