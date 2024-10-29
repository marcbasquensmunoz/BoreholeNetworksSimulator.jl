
function q_coef(::NoBoundary, medium, method, setup, i)
    constant_integral(medium, method, setup, i) + constant_coef(method, i)
end

function q_coef(::DirichletBoundaryCondition, medium, method, setup, i)
    @unpack expΔt, w, ζ = method
    constant_integral(medium, method, setup, i) - constant_integral(medium, method, image(setup), i) + constant_coef(method, i)
end

function q_coef(::NeumannBoundaryCondition, medium, method, setup, i)
    @unpack expΔt, w, ζ = method
    constant_integral(medium, method, setup, i) + constant_integral(medium, method, image(setup), i) + constant_coef(method, i)
end

function constant_coef(method::NonHistoryMethod, i)
    @unpack expΔt, w, ζ, aux = method
    @. aux = expΔt / ζ
    @views -dot(w[:, i], aux)
end

function constant_integral(medium::GroundMedium, method, setup::SegmentToSegment, i)
    @unpack D1, H1, D2, H2, σ = setup
    @unpack λ = medium

    β(d) = sqrt(σ^2 + d^2) + d*log(sqrt(σ^2 + d^2) - d)
    1/(4π*λ*H2) * (β(D1+H1-D2-H2) + β(D1-D2) - β(D1+H1-D2) - β(D1-D2-H2))
end

function constant_integral(medium::GroundMedium, method, setup::SegmentToPoint, i)
    @unpack D, H, z, σ = setup
    @unpack expΔt, w, ζ = method
    @unpack λ = medium

    1/(4π*λ) * log((z-D+sqrt(σ^2+(z-D)^2)) / (z-D-H+sqrt(σ^2+(z-D-H)^2)))
end
