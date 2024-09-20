using .FiniteLineSource: SegmentToSegment, SegmentToPoint, Constants, adaptive_gk_segments, DiscretizationParameters, f_guess, precompute_coefficients, IntegrationSegment

"""
    NonHistoryMethod{T} <: TimeSuperpositionMethod 
    NonHistoryMethod()

Use the non-history method to compute the thermal response between boreholes. 
See _A non-history dependent temporal superposition algorithm for the finite line source solution_.
It should be initialized without arguments, but it contains the variables:
- `F::Matrix{T}`: each column contains the `F` function (encoding the load history) for each borehole. It is initially 0.
- `ζ::Vector{T}`: discretization nodes of the integration interval. Shared for all boreholes. Precomputed in [`initialize`](@ref).
- `w::Matrix{T}`: weights of the ζ integration for each pair of boreholes. Precomputed in [`initialize`](@ref).
- `expΔt::Vector{T}`: exp(-ζ^2*Δt). Precomputed in [`initialize`](@ref).

This feature is experimental and might not work as expected in some cases. 
"""
mutable struct NonHistoryMethod{T} <: TimeSuperpositionMethod 
    F::Matrix{T}
    ζ::Vector{T}
    w::Matrix{T}
    expΔt::Vector{T}
    n_disc
    aux::Vector{T}
end
NonHistoryMethod(;n_disc=20) = NonHistoryMethod(zeros(0, 0), zeros(0), zeros(0, 0), zeros(0), n_disc, zeros(0))

FiniteLineSource.SegmentToSegment(s::MeanSegToSegEvParams) = SegmentToSegment(D1=s.D1, H1=s.H1, D2=s.D2, H2=s.H2, σ=s.σ)
image(s::SegmentToSegment) = SegmentToSegment(D1=-s.D1, H1=-s.H1, D2=s.D2, H2=s.H2, σ=s.σ)
image(s::SegmentToPoint) = SegmentToPoint(D=-s.D, H=-s.H, z=s.z, σ=s.σ)

function precompute_auxiliaries!(method::NonHistoryMethod, options)
    @unpack Nb, Nt, Ns, Δt, borefield, medium, boundary_condition = options
    @unpack n_disc = method
    α = get_α(medium)
    rb = get_rb(borefield, 1) 
    kg = get_λ(medium)
    Δt̃ = α*Δt/rb^2
    ϵ = 10^-18
    b = ceil(erfcinv(ϵ / sqrt(π/Δt̃)) / sqrt(Δt̃))

    constants = Constants(Δt=Δt, α=α, rb=rb, kg=kg, b=b)
    segments = [IntegrationSegment(0., 0.05, 0., 0.), IntegrationSegment(0.05, 0.1, 0., 0.), IntegrationSegment(0.1, 0.5, 0., 0.), IntegrationSegment(0.5, 1., 0., 0.), IntegrationSegment(1., 3., 0., 0.), IntegrationSegment(3., 7., 0., 0.)]
    #segments=adaptive_gk_segments(f_guess(SegmentToSegment(get_sts(borefield, 1, 1)), constants), 0., b)
    dps = [DiscretizationParameters(s.a, s.b, n_disc) for s in segments]
    ζ = reduce(vcat, (dp.x for dp in dps)) 
    expΔt = @. exp(-ζ^2 * Δt̃)

    n = length(ζ)
    w = zeros(n, Ns*Ns)

    containers, map = FiniteLineSource.initialize_containers(SegmentToSegment(get_sts(borefield, 1, 1)), dps)    

    for i in 1:Ns
        for j in 1:Ns
            setup = SegmentToSegment(get_sts(borefield, i, j))
            #setup = get_stp(borefield, i, j)
            w[:, (i-1)*Ns+j] = reduce(vcat, [coefficients(boundary_condition, setup, constants, dp, containers[map[k]]) for (k, dp) in enumerate(dps)])
        end
    end

    perm = sortperm(ζ)

    method.F = zeros(n, Ns*Ns)
    method.ζ = ζ[perm]
    method.w = w[perm, :]
    method.expΔt = expΔt[perm]
    method.aux = zeros(n)
end

function coefficients(::NoBoundary, setup::SegmentToSegment, params::Constants, dp, containers)
    precompute_coefficients(setup, params=params, dp=dp, containers=containers)
end

function coefficients(::DirichletBoundaryCondition, setup::SegmentToSegment, params::Constants, dp, containers)
    image_setup = image(setup)
    w1 = precompute_coefficients(setup, params=params, dp=dp, containers=containers)
    w2 = precompute_coefficients(image_setup, params=params, dp=dp, containers=containers)
    w1 - w2
end

function get_sts(borefield::Borefield, i, j)
    xi, yi, Di, Hi = segment_coordinates(borefield, i)
    xj, yj, Dj, Hj = segment_coordinates(borefield, j)
    σ = i == j ? get_rb(borefield, i) : sqrt((xi-xj)^2 + (yi-yj)^2)
    MeanSegToSegEvParams(D1=Di, H1=Hi, D2=Dj, H2=Hj, σ=σ)
end

function get_stp(borefield::Borefield, i, j)
    xi, yi, Di, Hi = segment_coordinates(borefield, i)
    xj, yj, Dj, Hj = segment_coordinates(borefield, j)
    σ = i == j ? get_rb(borefield, i) : sqrt((xi-xj)^2 + (yi-yj)^2)
    FiniteLineSource.SegmentToPoint(σ = σ, D = Di, H = Hi, z = Dj + Hj/2)
end

function update_auxiliaries!(method::NonHistoryMethod, X, borefield, step)
    @unpack ζ, F, expΔt = method
    Nb = n_boreholes(borefield)

    for i in 1:size(F)[2]
        @. @views F[:, i] = expΔt * F[:, i] + X[3Nb+(i-1)%Nb+1, step] * (1 - expΔt) / ζ
    end
end

function method_coeffs!(M, method::NonHistoryMethod, borefield, medium, boundary_condition)
    Nb = n_boreholes(borefield)
    Ns = n_segments(borefield)
    λ = get_λ(medium)

    for i in 1:Ns
        for j in 1:Ns
            setup = SegmentToSegment(get_sts(borefield, i, j))
            M[i, 3Nb+j] = q_coef(boundary_condition, medium, method, setup, λ, (i-1)*Ns+j) 
        end
    end

    for i in 1:Ns
        bh = where_is_segment(borefield, i)
        M[i, 2Nb + bh] = -1
    end
end

function method_b!(b, method::NonHistoryMethod, borefield, medium, step)
    @unpack w, expΔt, F, aux = method
    b .= -get_T0(medium)
    Nb = n_boreholes(borefield)

    @inbounds for i in eachindex(b)
        @inbounds for j in 1:Nb
            @views @. aux = expΔt * F[:, Nb*(i-1)+j]
            @views b[i] -= dot(w[:, Nb*(i-1)+j], aux)
        end
    end
end

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
