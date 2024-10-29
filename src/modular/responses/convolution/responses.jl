
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

function response(::NoBoundary, setup, params::Constants, t; atol, rtol)
    step_response(setup, params, t, atol=atol, rtol=rtol)
end

function response(::DirichletBoundaryCondition, setup, params::Constants, t; atol, rtol)
    setup_image = image(setup)
    Ip = @time step_response(setup, params, t, atol=atol, rtol=rtol)
    In = @time step_response(setup_image, params, t, atol=atol, rtol=rtol)
    Ip - In
end

function response(::AdiabaticBoundaryCondition, setup, params::Constants, t; atol, rtol)
    setup_image = image(setup)
    Ip = step_response(setup, params, t, atol=atol, rtol=rtol)
    In = step_response(setup_image, params, t, atol=atol, rtol=rtol)
    Ip + In
end

function response(::DirichletBoundaryCondition, s::SegmentToSegment, params::Constants, t; atol, rtol)
    @unpack rb, α, kg = params
    @unpack σ = s
    1 / (4π *s.H2 * kg) * quadgk(x -> exp(-σ^2*x^2) / x^2 * I_fls_dir(s, x), 1/sqrt(4*α*t), Inf, rtol = rtol, atol = atol)[1]
end

function response(::AdiabaticBoundaryCondition, s::SegmentToSegment, params::Constants, t; atol, rtol)
    @unpack rb, α, kg = params
    @unpack σ = s
    1 / (4π *s.H2 * kg) * quadgk(x -> exp(-σ^2*x^2) / x^2 * I_fls_neu(s, x), 1/sqrt(4*α*t), Inf, rtol = rtol, atol = atol)[1]
end

step_response(s::SegmentToPoint, params::Constants, t; atol, rtol) = stp_response(s, params, t, atol=atol, rtol=rtol)
step_response(s::SegmentToSegment, params::Constants, t; atol, rtol) = sts_response(s, params, t, atol=atol, rtol=rtol)

step_response(s::MovingSegmentToPoint, params::Constants, t; atol, rtol) = mstp_response(s, params, t, atol=atol, rtol=rtol)
step_response(s::MovingSegmentToSegment, params::Constants, t; atol, rtol) = msts_response(s, params, t, atol=atol, rtol=rtol)

function stp_response(s::SegmentToPoint, params::Constants, t; atol, rtol)
    @unpack σ, H, D, z = s
    @unpack α, kg = params
    1 / (4π * kg) * quadgk(x -> exp(-σ^2*x^2) / x * (erf(x * (z-D)) - erf(x * (z-D-H))), 1/sqrt(4*α*t), Inf, rtol = rtol, atol = atol)[1]
end

function sts_response(s::SegmentToSegment, params::Constants, t; atol, rtol)
    @unpack rb, α, kg = params
    @unpack σ = s
    1 / (4π *s.H2 * kg) * quadgk(x -> exp(-σ^2*x^2) / x^2 * I_fls(s, x), 1/sqrt(4*α*t), Inf, rtol = rtol, atol = atol)[1]
end

function mstp_response(s::MovingSegmentToPoint, params::Constants, t; atol, rtol)
    @unpack x, y, H, D, z, v = s
    @unpack α, kg = params
    return mfls_step_response(t, x, y, z, v, D, D+H, α, kg, atol=atol, rtol=rtol)
end

function msts_response(s::MovingSegmentToSegment, params::Constants, t; atol, rtol)
    @unpack α, kg = params
    @unpack x, v, D1, H1, D2, H2 = s
    params = FiniteLineSource.MeanSegToSegEvParams(s)
    r_min, r_max = FiniteLineSource.h_mean_lims(params)
    f(r) = FiniteLineSource.h_mean_sts(r, params) * moving_point_step_response(t, r, x, v, α, kg)
    X, W = FiniteLineSource.adaptive_nodes_and_weights(f, r_min, r_max, atol=atol, rtol=rtol)
    return dot(f.(X), W)
end

ierf(x) = x*erf(x) - (1- exp(-x^2)) / sqrt(π)
I_fls(s::SegmentToSegment, x) = ierf((s.D2 - s.D1 + s.H2)*x) + ierf((s.D2 - s.D1 - s.H1)*x) - ierf((s.D2 - s.D1)*x) - ierf((s.D2 - s.D1 + s.H2 - s.H1)*x)
I_fls_dir(s::SegmentToSegment, x) = ierf((s.D2 - s.D1 + s.H2)*x) + ierf((s.D2 - s.D1 - s.H1)*x) - ierf((s.D2 - s.D1)*x) - ierf((s.D2 - s.D1 + s.H2 - s.H1)*x) + ierf((s.D2 + s.D1 + s.H2)*x) + ierf((s.D2 + s.D1 + s.H1)*x) - ierf((s.D2 + s.D1)*x) - ierf((s.D2 + s.D1 + s.H2 + s.H1)*x)
I_fls_neu(s::SegmentToSegment, x) = ierf((s.D2 - s.D1 + s.H2)*x) + ierf((s.D2 - s.D1 - s.H1)*x) - ierf((s.D2 - s.D1)*x) - ierf((s.D2 - s.D1 + s.H2 - s.H1)*x) - ierf((s.D2 + s.D1 + s.H2)*x) - ierf((s.D2 + s.D1 + s.H1)*x) + ierf((s.D2 + s.D1)*x) + ierf((s.D2 + s.D1 + s.H2 + s.H1)*x)

point_step_response(t, r, α, kg) = erfc(r/(2*sqrt(t*α))) / (4*π*r*kg)   
moving_point_step_response(t, r, x, v, α, kg) = exp(-v * (r-x)/ (2α)) * (erfc( (r-t*v) / sqrt(4t*α)) + erfc((r+t*v) / sqrt(4t*α)) * exp(v*r/α) ) / (8π*r*kg)

mfls_step_response(t, x, y, z, v, a, b, α, kg; atol, rtol) = quadgk(ζ -> moving_point_step_response(t, sqrt(x^2 + y^2 + (z-ζ)^2), x, v, α, kg), a, b, atol=atol, rtol=rtol)[1]
