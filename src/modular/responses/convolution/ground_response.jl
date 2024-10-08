
function compute_response!(g, medium::GroundMedium, borefield, boundary_condition, t) 
    @unpack λ, α = medium
    Ns = n_segments(borefield)

    coords = [segment_coordinates(borefield, i) for i in 1:Ns]
    
    for (k, tt) in enumerate(t)
        for j in 1:Ns
            for i in 1:Ns
                x1, y1, D1, H1 = coords[i]
                x2, y2, D2, H2 = coords[j]
                rb = get_rb(borefield, i)
                σ = i == j ? rb : sqrt((x2-x1)^2 + (y2-y1)^2)
                s = SegmentToSegment(D1=D1, H1=H1, D2=D2, H2=H2, σ=σ)
                g[i, j, k] = sts(boundary_condition, s, Constants(α=α, kg=λ, rb=rb), tt)
            end
        end
    end
end

point_step_response(t, r, α, kg) = erfc(r/(2*sqrt(t*α))) / (4*π*r*kg)    
fls_step_response(t, σ, z, a, b, α, kg; atol=1e-8) = quadgk(ζ -> point_step_response(t, sqrt(σ^2 + (z-ζ)^2), α, kg), a, b, atol=atol)[1]

function fls_adiabatic_surface(t, x, y, z, H, D, α, kg, atol = 1e-8)
    σ = sqrt(x^2+y^2)
    I1 = fls_step_response(t, σ, z, D, D+H, α, kg, atol=atol)
    I2 = fls_step_response(t, σ, z, -D-H, -D, α, kg, atol=atol) 
    return (I1 + I2)
end

function fls_dirichlet(t, x, y, z, H, D, α, kg, atol = 1e-8)
    σ = sqrt(x^2+y^2)
    I1 = fls_step_response(t, σ, z, D, D+H, α, kg, atol=atol)
    I2 = fls_step_response(t, σ, z, -D-H, -D, α, kg, atol=atol) 
    return (I1 - I2)
end

function fls(t, x, y, z, H, D, α, kg, atol = 1e-8)
    σ = sqrt(x^2+y^2)
    return fls_step_response(t, σ, z, D, D+H, α, kg, atol=atol)
end

function sts_response(s::SegmentToSegment, params::Constants, t)
    @unpack α, kg = params
    params = FiniteLineSource.MeanSegToSegEvParams(s)
    r_min, r_max = FiniteLineSource.h_mean_lims(params)
    f(r) = FiniteLineSource.h_mean_sts(r, params) * point_step_response(t, r, α, kg)
    x, w = FiniteLineSource.adaptive_nodes_and_weights(f, r_min, r_max)
    dot(f.(x), w)
end

function sts(::NoBoundary, s::SegmentToSegment, params::Constants, t)
    sts_response(s, params, t)
end

function sts(::DirichletBoundaryCondition, s::SegmentToSegment, params::Constants, t)
    s_image = SegmentToSegment(D1=-s.D1, H1=-s.H1, D2=s.D2, H2=s.H2, σ=s.σ)
    Ip = sts_response(s, params, t)
    In = sts_response(s_image, params, t)
    Ip - In
end
