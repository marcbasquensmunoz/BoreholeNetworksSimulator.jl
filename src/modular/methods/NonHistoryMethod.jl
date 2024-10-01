using .FiniteLineSource: SegmentToSegment, SegmentToPoint, Constants, DiscretizationParameters, f_guess, precompute_coefficients

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
    n_disc::Int
    aux::Vector{T}
end
NonHistoryMethod(;n_disc::Int=20) = NonHistoryMethod(zeros(0, 0), zeros(0), zeros(0, 0), zeros(0), n_disc, zeros(0))

function precompute_auxiliaries!(method::NonHistoryMethod, options)
    @unpack Nb, Nt, Ns, Δt, borefield, medium, boundary_condition, approximation = options
    @unpack n_disc = method
    α = get_α(medium)
    rb = get_rb(borefield, 1) 
    kg = get_λ(medium)
    Δt̃ = α*Δt/rb^2
    ϵ = 10^-18
    b = ceil(erfcinv(ϵ / sqrt(π/Δt̃)) / sqrt(Δt̃))

    constants = Constants(Δt=Δt, α=α, rb=rb, kg=kg, b=b)
    _, _, segments = quadgk_segbuf(f_guess(setup(approximation, borefield, 1, 1), constants), 0., b)
    @show length(segments)
    dps = [DiscretizationParameters(s.a, s.b, n_disc) for s in segments]
    ζ = reduce(vcat, (dp.x for dp in dps)) 
    expΔt = @. exp(-ζ^2 * Δt̃)

    n = length(ζ)
    w = zeros(n, Ns*Ns)

    containers, map = FiniteLineSource.initialize_containers(setup(approximation, borefield, 1, 1), dps)    

    for i in 1:Ns
        for j in 1:Ns
            w[:, (i-1)*Ns+j] = reduce(vcat, [weights(boundary_condition, setup(approximation, borefield, i, j), constants, dp, containers[map[k]]) for (k, dp) in enumerate(dps)])
        end
    end

    perm = sortperm(ζ)

    method.F = zeros(n, Ns*Ns)
    method.ζ = ζ[perm]
    method.w = w[perm, :]
    method.expΔt = expΔt[perm]
    method.aux = zeros(n)
end

function update_auxiliaries!(method::NonHistoryMethod, X, borefield, step)
    @unpack ζ, F, expΔt = method
    Nb = n_boreholes(borefield)

    for i in 1:size(F)[2]
        @. @views F[:, i] = expΔt * F[:, i] + X[3Nb+(i-1)%Nb+1, step] * (1 - expΔt) / ζ
    end
end

function method_coeffs!(M, method::NonHistoryMethod, borefield, medium, boundary_condition, approximation)
    Nb = n_boreholes(borefield)
    Ns = n_segments(borefield)
    λ = get_λ(medium)

    for i in 1:Ns
        for j in 1:Ns
            M[i, 3Nb+j] = q_coef(boundary_condition, medium, method, setup(approximation, borefield, i, j), λ, (i-1)*Ns+j) 
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
