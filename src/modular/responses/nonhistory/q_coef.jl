
function q_coef(::NoBoundary, medium, method, setup, i)
    (method isa NonHistoryMethod && (i-1)%2+1 != div(i-1, 2)+1 ? 0. : constant_integral(medium, setup, i)) + constant_coef(method, i)
end

function q_coef(::DirichletBoundaryCondition, medium, method, setup, i)
    constant_integral(medium, setup, i) - constant_integral(medium, image(setup), i) + constant_coef(method, i)
end

function q_coef(::NeumannBoundaryCondition, medium, method, setup, i)
    constant_integral(medium, setup, i) + constant_integral(medium, image(setup), i) + constant_coef(method, i)
end

function constant_coef(method::OriginalNonHistoryMethod, i)
    @unpack expΔt, w, ζ, aux = method
    @. aux = expΔt / ζ
    @show -dot(w[:, i], aux)
    @views -dot(w[:, i], aux)
end

function constant_coef(method::NonHistoryMethod, i)
    @unpack HM, sr_F, sr_w, sr_expt, sr_ζ = method
    Nb = size(HM, 2)
    source = (i-1)%Nb+1
    target = div(i-1, Nb)+1

    res = 0.
    if source == target
        res = -dot(sr_expt[source] ./ sr_ζ[source], sr_w[source])
    end
    return res
end

function constant_integral(medium::GroundMedium, setup::SegmentToSegment, i)
    @unpack D1, H1, D2, H2, σ = setup
    @unpack λ = medium

    β(d) = sqrt(σ^2 + d^2) + d*log(sqrt(σ^2 + d^2) - d)
    1/(4π*λ*H2) * (β(D1+H1-D2-H2) + β(D1-D2) - β(D1+H1-D2) - β(D1-D2-H2))
end

function constant_integral(medium::GroundMedium, setup::SegmentToPoint, i)
    @unpack D, H, z, σ = setup
    @unpack λ = medium

    1/(4π*λ) * log((z-D+sqrt(σ^2+(z-D)^2)) / (z-D-H+sqrt(σ^2+(z-D-H)^2)))
end
