using Parameters

@with_kw struct GroundMedium <: Medium
    λ                
    C
    α = λ/C 
end
get_λ(bfp::GroundMedium) = bfp.λ

function compute_response!(medium::GroundMedium, borefield::Borefield, coord_source, coord_eval, t) end