
function compute_response!(g, medium::Medium, options) 
    @unpack λ, α = medium
    @unpack borefield, boundary_condition, approximation, atol, rtol, t = options
    Nb = n_boreholes(borefield)

    for (k, tt) in enumerate(t)
        for j in 1:Nb
            for i in 1:Nb
                rb = get_rb(borefield, i)
                s = setup(approximation, medium, borefield, i, j)
                g[i, j, k] = response(boundary_condition, s, Constants(α=α, kg=λ, rb=rb), tt, atol=atol, rtol=rtol)
            end
        end
    end
end

# TODO: Might be able to improve performance by explicitly writing the integrand in each case
integrand(::NoBoundary, setup, s) = integrand(setup, s)
integrand(::DirichletBoundaryCondition, setup, s) = integrand(setup, s) - integrand(image(setup), s)
integrand(::NeumannBoundaryCondition, setup, s) = integrand(setup, s) + integrand(image(setup), s)

integrand(setup::SegmentToPoint, s) = I_stp(s, setup.D, setup.H, setup.z)
integrand(setup::SegmentToSegment, s) = I_sts(s, setup.D1, setup.D2, setup.H1, setup.H2) 
integrand(setup::MovingSegmentToPoint, s) = I_stp(s, setup.D, setup.H, setup.z)
integrand(setup::MovingSegmentToSegment, s) = I_sts(s, setup.D1, setup.D2, setup.H1, setup.H2) 

function response(boundary_condition, setup::SegmentToPoint, params::Constants, t; atol, rtol)
    @unpack α, kg = params
    @unpack σ, H, D, z = setup
    1 / (4π * kg) * quadgk(s -> exp(-σ^2*s^2) / s * integrand(boundary_condition, setup, s), 1/sqrt(4*α*t), Inf, rtol = rtol, atol = atol)[1]
end

function response(boundary_condition, setup::SegmentToSegment, params::Constants, t; atol, rtol)
    @unpack α, kg = params
    @unpack σ, H2 = setup
    1 / (4π * H2 * kg) * quadgk(s -> exp(-σ^2*s^2) / s^2 * integrand(boundary_condition, setup, s), 1/sqrt(4*α*t), Inf, rtol = rtol, atol = atol)[1]
end

function response(boundary_condition, setup::MovingSegmentToPoint, params::Constants, t; atol, rtol)
    @unpack α, kg = params
    @unpack x, y, H, D, z, v = setup
    1 / (4π * kg) * quadgk(s -> exp(-((x-v/(4α*s^2))^2 + y^2)*s^2) / s * integrand(boundary_condition, setup, s), 1/sqrt(4*α*t), Inf, rtol = rtol, atol = atol)[1]
end

function response(boundary_condition, setup::MovingSegmentToSegment, params::Constants, t; atol, rtol)
    @unpack α, kg = params
    @unpack x, y, v, D1, H1, D2, H2 = setup
    1 / (4π * H2 * kg) * quadgk(s -> exp(-s^2*((x - v/(4α*s^2))^2 + y^2)) / s^2 * integrand(boundary_condition, setup, s), 1/sqrt(4*α*t), Inf, rtol = rtol, atol = atol)[1]
end

ierf(x) = x*erf(x) - (1 - exp(-x^2)) / sqrt(π)
I_sts(s, D1, D2, H1, H2) = ierf((D2 - D1 + H2)*s) + ierf((D2 - D1 - H1)*s) - ierf((D2 - D1)*s) - ierf((D2 - D1 + H2 - H1)*s)
I_stp(s, D, H, z) = erf(s * (z-D)) - erf(s * (z-D-H))
