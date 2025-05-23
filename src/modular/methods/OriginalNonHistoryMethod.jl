using .FiniteLineSource: SegmentToSegment, SegmentToPoint, Constants, DiscretizationParameters, f_guess, precompute_coefficients, initialize_containers, make_DiscretizationParameters

"""
    OriginalNonHistoryMethod{T} <: TimeSuperpositionMethod 
    OriginalNonHistoryMethod()

Use the non-history method to compute the thermal response between boreholes. 
See _A non-history dependent temporal superposition algorithm for the finite line source solution_.
It should be initialized without arguments, but it contains the variables:
- `F::Matrix{T}`: each column contains the `F` function (encoding the load history) for each borehole. It is initially 0.
- `ζ::Vector{T}`: discretization nodes of the integration interval. Shared for all boreholes. Precomputed in [`initialize`](@ref).
- `w::Matrix{T}`: weights of the ζ integration for each pair of boreholes. Precomputed in [`initialize`](@ref).
- `expΔt::Vector{T}`: exp(-ζ^2*Δt). Precomputed in [`initialize`](@ref).

This feature is experimental and might not work as expected in some cases. 
"""
mutable struct OriginalNonHistoryMethod{T} <: TimeSuperpositionMethod 
    F::Matrix{T}
    ζ::Vector{T}
    w::Matrix{T}
    expΔt::Vector{T}
    n_disc::Int
    aux::Vector{T}
end
OriginalNonHistoryMethod(;n_disc::Int=20) = OriginalNonHistoryMethod(zeros(0, 0), zeros(0), zeros(0, 0), zeros(0), n_disc, zeros(0))

function distances(borefield, boundary_condition, approximation, medium)
    map = Dict{setup_type(approximation, medium), Int}()
    buffers = get_buffers(boundary_condition)

    k = 1
    boreholes = 1:n_boreholes(borefield)
    for i in boreholes
        for j in boreholes
            s = setup(approximation, medium, borefield, i, j)
            if !haskey(map, s)
                map[s] = k
                k += 1

                rb = get_rb(borefield, i)
                add_buffer!(buffers, boundary_condition, s, rb)
            end
        end
    end
    return map, buffers
end

function precompute_auxiliaries!(method::OriginalNonHistoryMethod, options)
    @unpack Nb, Nt, Ns, Δt, borefield, medium, boundary_condition, approximation, atol, rtol = options
    @unpack n_disc = method
    α = get_α(medium)
    rb = get_rb(borefield, 1) 
    kg = get_λ(medium)
    Δt̃ = α*Δt/rb^2
    ϵ = 10^-18
    b = ceil(erfcinv(ϵ / sqrt(π/Δt̃)) / sqrt(Δt̃))

    constants = Constants(Δt=Δt, α=α, rb=rb, kg=kg, b=b)
    _, _, segments = quadgk_segbuf(f_guess(setup(approximation, medium, borefield, 1, 1), constants), 0., b, atol=atol, rtol=rtol)
    xt, wt = gausslegendre(n_disc+1)  
    dps = [make_DiscretizationParameters(s.a, s.b, n_disc, xt=xt, w=wt) for s in segments]
    ζ = reduce(vcat, (dp.x for dp in dps)) 
    expΔt = @. exp(-ζ^2 * Δt̃)

    distances_map, quadgk_buffers = distances(borefield, boundary_condition, approximation, medium)
    disc_map, containers = initialize_containers(setup(approximation, medium, borefield, 1, 1), dps)    

    n = length(ζ)
    w = zeros(n, Ns*Ns)
    w_buffer = zeros(n, length(distances_map))

    for (key, value) in pairs(distances_map)
        for (k, dp) in enumerate(dps)
            range = (n_disc+1)*(k-1)+1:(n_disc+1)*k
            w_buffer[range, value] .= weights(boundary_condition, key, constants, dp, containers[disc_map[k]], quadgk_buffers[value], atol=atol, rtol=rtol)
        end
    end
    for i in 1:Ns
        for j in 1:Ns
            k = distances_map[setup(approximation, medium, borefield, i, j)]
            @inbounds @views @. w[:, (i-1)*Ns+j] = w_buffer[:, k]
        end
    end

    perm = sortperm(ζ)

    @views method.ζ = ζ[perm]
    @views method.w = w[perm, :]
    @views method.expΔt = expΔt[perm]
    method.F = zeros(n, Ns*Ns)
    method.aux = zeros(n)
end

function update_auxiliaries!(method::OriginalNonHistoryMethod, X, borefield, step)
    @unpack ζ, F, expΔt = method
    Nb = n_boreholes(borefield)

    for i in 1:size(F)[2]
        @. @views F[:, i] = expΔt * F[:, i] + X[3Nb+(i-1)%Nb+1, step] * (1 - expΔt) / ζ
    end
end

function method_coeffs!(M, method::OriginalNonHistoryMethod, options)
    @unpack borefield, medium, boundary_condition, approximation = options
    Nb = n_boreholes(borefield)
    Ns = n_segments(borefield)

    for i in 1:Ns
        for j in 1:Ns
            M[i, 3Nb+j] = q_coef(boundary_condition, medium, method, setup(approximation, medium, borefield, i, j), (i-1)*Ns+j) 
        end
    end

    for i in 1:Ns
        bh = where_is_segment(borefield, i)
        M[i, 2Nb + bh] = -1
    end
end

function method_b!(b, method::OriginalNonHistoryMethod, borefield, medium, step)
    @unpack w, expΔt, F, aux = method
    b .= -get_T0(medium)
    Nb = n_boreholes(borefield)

    @inbounds for i in eachindex(b)
        @inbounds for j in 1:Nb
            @views @. aux = expΔt * F[:, Nb*(i-1)+j]
            @show i, j, dot(w[:, Nb*(i-1)+j], aux)
            @views b[i] -= dot(w[:, Nb*(i-1)+j], aux)
        end
    end
end

function reset!(method::OriginalNonHistoryMethod) 
    method.F .= 0
    method.aux .= 0
end
