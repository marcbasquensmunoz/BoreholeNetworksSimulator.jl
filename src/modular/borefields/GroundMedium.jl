using Parameters

@with_kw struct GroundMedium <: Medium
    λ                
    C
    α = λ/C 
end
λ(bfp::GroundMedium) = bfp.λ
