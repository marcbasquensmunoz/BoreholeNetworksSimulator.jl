using .FiniteLineSource: SegmentToSegment, SegmentToPoint, Constants, DiscretizationParameters, f_guess, precompute_coefficients, initialize_containers

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

function get_boreholes_distance(borefield, i, j)
    x1, y1, D1, H1 = segment_coordinates(borefield, i)
    x2, y2, D2, H2 = segment_coordinates(borefield, j)

    σ = i == j ? get_rb(borefield, i) : sqrt((x1 - x2)^2 + (y1 - y2)^2)
    return SegmentToSegment(σ=σ, D1=D1, D2=D2, H1=H1, H2=H2)
end

function distances(borefield, boundary_condition)
    map = Dict{SegmentToSegment{Float64}, Int}()
    buffers = get_buffers(boundary_condition)

    k = 1
    boreholes = 1:n_boreholes(borefield)
    for i in boreholes
        for j in boreholes
            s = get_boreholes_distance(borefield, i, j)
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
    dps = [DiscretizationParameters(s.a, s.b, n_disc) for s in segments]
    ζ = reduce(vcat, (dp.x for dp in dps)) 
    expΔt = @. exp(-ζ^2 * Δt̃)

    distances_map, quadgk_buffers = distances(borefield, boundary_condition)
    disc_map, containers = initialize_containers(setup(approximation, borefield, 1, 1), dps)    

    n = length(ζ)
    w = zeros(n, Ns*Ns)
    w_buffer = zeros(n, length(distances_map))

    for (key, value) in pairs(distances_map)
        for (k, dp) in enumerate(dps)
            range = (n_disc+1)*(k-1)+1:(n_disc+1)*k
            w_buffer[range, value] .= weights(boundary_condition, key, constants, dp, containers[disc_map[k]], quadgk_buffers[value])
        end
    end
    for i in 1:Ns
        for j in 1:Ns
            k = distances_map[get_boreholes_distance(borefield, i, j)]
            @views @. w[:, (i-1)*Ns+j] = w_buffer[:, k]
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

function method_coeffs!(M, method::NonHistoryMethod, options)
    @unpack borefield, medium, boundary_condition, approximation = options
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
