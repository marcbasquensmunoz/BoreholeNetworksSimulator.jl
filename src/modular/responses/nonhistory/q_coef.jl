
function q_coef(::NoBoundary, medium, method, setup, λ, i)
    constant_integral(medium, method, setup, λ, i) + constant_coef(method, i)
end

function q_coef(::DirichletBoundaryCondition, medium, method, setup, λ, i)
    @unpack expΔt, w, ζ = method
    constant_integral(medium, method, setup, λ, i) - constant_integral(medium, method, image(setup), λ, i) + constant_coef(method, i)
end

function constant_coef(method::NonHistoryMethod, i)
    @unpack expΔt, w, ζ, aux = method
    @. aux = expΔt / ζ
    @views -dot(w[:, i], aux)
end

function constant_integral(::GroundMedium, method, setup::SegmentToSegment, λ, i)
    @unpack D1, H1, D2, H2, σ = setup

    β(d) = sqrt(σ^2 + d^2) + d*log(sqrt(σ^2 + d^2) - d)
    1/(4π*λ*H2) * (β(D1+H1-D2-H2) + β(D1-D2) - β(D1+H1-D2) - β(D1-D2-H2))
end

function constant_integral(::GroundMedium, method, setup::SegmentToPoint, λ, i)
    @unpack D, H, z, σ = setup
    @unpack expΔt, w, ζ = method

    1/(4π*λ) * log((z-D+sqrt(σ^2+(z-D)^2))/(z-D-H+sqrt(σ^2+(z-D-H)^2)))
end
