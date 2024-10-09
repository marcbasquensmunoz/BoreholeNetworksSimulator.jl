
function compute_response!(g, medium::Medium, borefield, boundary_condition, approximation, t) 
    @unpack λ, α = medium
    Nb = n_boreholes(borefield)

    for (k, tt) in enumerate(t)
        for j in 1:Nb
            for i in 1:Nb
                rb = get_rb(borefield, i)
                s = setup(approximation, borefield, i, j)
                g[i, j, k] = response(boundary_condition, medium, s, Constants(α=α, kg=λ, rb=rb), tt)
            end
        end
    end
end

function response(::NoBoundary, medium, setup, params::Constants, t)
    step_response(setup, medium, params, t)
end

function response(::DirichletBoundaryCondition, medium, setup, params::Constants, t)
    setup_image = image(setup)
    Ip = step_response(setup, medium, params, t)
    In = step_response(setup_image, medium, params, t)
    Ip - In
end

step_response(s::SegmentToPoint, ::GroundMedium, params::Constants, t) = stp_response(s, params, t)
step_response(s::SegmentToSegment, ::GroundMedium, params::Constants, t) = sts_response(s, params, t)

# Fix this
step_response(s::SegmentToPoint, ::FlowInPorousMedium, params::Constants, t) = stp_response(s, params, t)
step_response(s::SegmentToSegment, ::FlowInPorousMedium, params::Constants, t) = sts_response(s, params, t)

function sts_response(s::SegmentToSegment, params::Constants, t)
    @unpack α, kg = params
    params = FiniteLineSource.MeanSegToSegEvParams(s)
    r_min, r_max = FiniteLineSource.h_mean_lims(params)
    f(r) = FiniteLineSource.h_mean_sts(r, params) * point_step_response(t, r, α, kg)
    x, w = FiniteLineSource.adaptive_nodes_and_weights(f, r_min, r_max)
    return dot(f.(x), w)
end

function stp_response(s::SegmentToPoint, params::Constants, t)
    @unpack σ, H, D, z = s
    @unpack α, kg = params
    return fls_step_response(t, σ, z, D, D+H, α, kg)
end

fls_step_response(t, σ, z, a, b, α, kg; atol=1e-8) = quadgk(ζ -> point_step_response(t, sqrt(σ^2 + (z-ζ)^2), α, kg), a, b, atol=atol)[1]
point_step_response(t, r, α, kg) = erfc(r/(2*sqrt(t*α))) / (4*π*r*kg)    
mfls_step_response(t, σ, z, a, b, α, kg) = 1/(2π*kg)*mfls_adiabatic_surface(t, α, x, y, z, vt, h, D; atol =1e-9) 
