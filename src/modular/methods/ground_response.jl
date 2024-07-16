
function compute_response!(g, medium::GroundMedium, borefield::Borefield, t) 
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
                g[i, j, k] = sts(s, Constants(α=α, kg=λ, rb=rb), tt)
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

function sts(s:: SegmentToSegment, params::Constants, t)
    @unpack α, kg = params
    params = FiniteLineSource.MeanSegToSegEvParams(s)
    h_mean_sts, r_min, r_max = FiniteLineSource.mean_sts_evaluation(params)
    f(r) = h_mean_sts(r) * point_step_response(t, r, α, kg)
    x, w = FiniteLineSource.adaptive_gk(f, r_min, r_max)
    dot(f.(x), w)
end
