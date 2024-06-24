
@with_kw struct GroundMedium{T <: Real} <: Medium @deftype T
    λ = 3.               
    α = 1e-6
    C = λ/α
end
get_λ(medium::GroundMedium) = medium.λ
get_α(medium::GroundMedium) = medium.α

function compute_response!(g, medium::GroundMedium, borefield::Borefield, t) 
end