mutable struct NonHistoryMethod{T} <: Method 
    F::Matrix{T}
    ζ::Vector{T}
    w::Matrix{T}
    expΔt::Vector{T}
end
check_compatibility(borefield::Borefield, ::Constraint, ::NonHistoryMethod) = borefield.medium isa GroundWaterMedium ? NotCompatible("The non-history method is not implemented with ground water flow yet") : Compatible() 

FiniteLineSource.SegmentToSegment(s::MeanSegToSegEvParams) = FiniteLineSource.SegmentToSegment(D1=s.D1, H1=s.H1, D2=s.D2, H2=s.H2, σ=s.σ)

function NonHistoryMethod(;parameters, borefield, b = 10.)
    @unpack Nb, Nt, Ns, tstep = parameters
    Δt = tstep
    α = get_α(borefield.medium)
    rb = get_rb(borefield, 1) 
    kg = get_λ(borefield.medium)
    Δt̃ = α*Δt/rb^2

    n_segment = 20

    constants = FiniteLineSource.Constants(Δt=Δt, α=α, rb=rb, kg=kg, b=b)
    segments = FiniteLineSource.adaptive_gk_segments(FiniteLineSource.f_guess(FiniteLineSource.SegmentToSegment(get_sts(borefield, 1, 1)), constants), 0., b)
    dps = @views [FiniteLineSource.discretization_parameters(s.a, s.b, n_segment) for s in segments]
    ζ = reduce(vcat, (dp.x for dp in dps)) 
    expΔt = @. exp(-ζ^2 * Δt̃)

    n = length(ζ)
    w = zeros(n, Ns*Ns)

    for i in 1:Ns
        for j in 1:Ns
            rb = get_rb(borefield, (i-1)*Ns+j) # Get the right rb
            sts = get_sts(borefield, i, j)
            w[:, (i-1)*Ns+j] = reduce(vcat, [FiniteLineSource.precompute_coefficients(FiniteLineSource.SegmentToSegment(sts), params=constants, dp=dp) for dp in dps])
        end
    end

    perm = sortperm(ζ)

    return NonHistoryMethod(zeros(n, Ns), ζ[perm], w[perm, :], expΔt[perm])
end

function get_sts(borefield::Borefield, i, j)
    xi, yi, Di, Hi = segment_coordinates(borefield, i)
    xj, yj, Dj, Hj = segment_coordinates(borefield, j)
    σ = i == j ? get_rb(borefield, i) : sqrt((xi-xj)^2 + (yi-yj)^2)
    MeanSegToSegEvParams(D1=Di, H1=Hi, D2=Dj, H2=Hj, σ=σ)
end

function update_auxiliaries!(method::NonHistoryMethod, X, current_Q, borefield::Borefield, step)
    @unpack ζ, F, expΔt = method

    for i in eachindex(size(F)[2])
        @. @views F[:, i] = expΔt * F[:, i] + current_Q[step] * (1 - expΔt) / ζ
    end
end

function method_coeffs!(M, method::NonHistoryMethod, borefield::Borefield)
    Nb = borehole_amount(borefield)
    Ns = segment_amount(borefield)
    λ = get_λ(borefield.medium)

    for i in 1:Ns
        for j in 1:Ns
            sts = get_sts(borefield, i, j)
            rb = get_rb(borefield, i) 
            M[i, 3Nb+j] = - q_coef(borefield.medium, method, sts, λ, (i-1)*Ns+j)
        end
    end

    for i in 1:Ns
        bh = where_is_segment(borefield, i)
        M[i, 2Nb + bh] = 1
    end
end

function method_b!(b, method::NonHistoryMethod, borefield::Borefield, step, current_Q)
    @unpack w, expΔt, F = method
    λ = get_λ(borefield.medium)
    C = 1 / (2*λ^2*π^2)
    Nb = borehole_amount(borefield)
    b .= get_T0(borefield)

    @show current_Q

    for i in eachindex(b)
        sts = get_sts(borefield, div(i-1,Nb)+1, (i-1)%Nb+1)
        rb = get_rb(borefield, i) 
        @inbounds b[i] += C * dot(w[:, i], expΔt .* F[:, i]) + current_Q[i] * q_coef(borefield.medium, method, sts, λ, i)
    end
end

function q_coef(::GroundMedium, method, sts, λ, i)
    @unpack D1, H1, D2, H2, σ = sts
    @unpack expΔt, w, ζ = method

    β(d) = sqrt(σ^2 + d^2) + d*log(sqrt(σ^2 + d^2) - d)
    constant_integral = 1/(4π*λ*H2) * (β(D1+H1-D2-H2) + β(D1-D2) - β(D1+H1-D2) - β(D1-D2-H2))

    constant_integral - dot(w[:, i], expΔt ./ ζ)
end
