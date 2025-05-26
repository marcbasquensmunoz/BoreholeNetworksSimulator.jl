using FiniteLineSource: prepare_containers, LineSource, AsymptoticContainers

mutable struct NonHistoryMethod{T <: Number} <: TimeSuperpositionMethod 
    ζ::Vector{T}
    F::Matrix{T}
    expt::Vector{T}
    expNin::Vector{T}
    expNout::Vector{T}
    HM::Array{T, 3}
    ranges::Vector{UnitRange{Int}}
    Kranges::Vector{UnitRange{Int}}
    K_min::Matrix{Int}
    qaux::Vector{T}
    N::Vector{Int}
    sr_ζ::Vector{Vector{T}}
    sr_w::Vector{Vector{T}}
    sr_F::Vector{Vector{T}}
    sr_expt::Vector{Vector{T}}
    sr_expNout::Vector{Vector{T}}
    sr_Ic::Vector{T}
    sr_Icout::Vector{T}
    q_N1::Vector{T}
end
NonHistoryMethod() = NonHistoryMethod(zeros(0), zeros(0, 0), zeros(0), zeros(0), zeros(0), zeros(0, 0, 0), UnitRange{Int}[], UnitRange{Int}[], zeros(Int, 0, 0), zeros(0), zeros(Int, 0), Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], zeros(0), zeros(0), zeros(0))

get_setup(::MidPointApproximation) = SegmentToPoint(H=0., D=0., z=0., σ=0.)
get_setup(::MeanApproximation) = SegmentToSegment(H1=0., D1=0., H2=0., D2=0., σ=0.)

function precompute_auxiliaries!(method::NonHistoryMethod, options)
    constants = Constants(Δt=options.Δt, α=get_α(options.medium), rb=get_rb(options.borefield, 1), kg=get_λ(options.medium))
    containers = AsymptoticContainers(10)
    ϵ = options.atol != 0 ? options.atol : 1e-4

    Nb = n_boreholes(options.borefield)
    sources = [LineSource(segment_coordinates(options.borefield, i)..., get_rb(options.borefield, i)) for i in 1:Nb]

    rep_setup = FiniteLineSource.get_representative_ltl(sources)[1] 
    blocks = prepare_containers(rep_setup, sources, ϵ, options.Nt, constants, containers)

    method.F = blocks.F
    method.ζ = blocks.ζ
    method.expt = blocks.expt
    method.expNin = blocks.expNin
    method.expNout = blocks.expNout
    method.HM = blocks.HM
    method.ranges = blocks.ranges
    method.Kranges = blocks.Kranges
    method.K_min = blocks.K_min
    method.qaux = blocks.qaux
    method.N = blocks.N
    method.sr_ζ = blocks.sr_ζ
    method.sr_w = blocks.sr_w
    method.sr_F = blocks.sr_F
    method.sr_expt = blocks.sr_expt
    method.sr_expNout = blocks.sr_expNout
    method.sr_Ic = blocks.sr_Ic
    method.sr_Icout = blocks.sr_Icout
    method.q_N1 = zeros(Nb)
end

function update_auxiliaries!(method::NonHistoryMethod, X, borefield, step)
    @unpack sr_F, sr_expt, sr_expNout, sr_ζ, F, N, expt, expNin, expNout, qaux, ranges, q_N1 = method

    Nb = n_boreholes(borefield)
    K = length(N) - 1
    @views q = X[3Nb+1:4Nb, 1:step]

    for i in 1:Nb
        @. sr_F[i] = sr_expt[i] * sr_F[i] + (1 - sr_expt[i]) / sr_ζ[i] * (q[i, step] - sr_expNout[i] * q_N1[i])
    end
    @views @. q_N1 = step - N[1] + 1 > 0 ? q[:, step - N[1] + 1] : zeros(Nb)

    for j in 1:Nb
        qaux .= 0.
        for i in 1:K
            qin = step - N[i] + 1 > 0 ? q[j, step - N[i] + 1] : 0.
            qout = step - N[i+1] + 1 > 0 ? q[j, step - N[i+1] + 1] : 0.
            @views @inbounds @. qaux[ranges[i]] = qin * expNin[ranges[i]] - qout * expNout[ranges[i]]
        end
        @views @. F[:, j] = expt * F[:, j] + qaux
    end
end

function method_coeffs!(M, method::NonHistoryMethod, options)
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

function method_b!(b, method::NonHistoryMethod, borefield, medium, step)
    @unpack F, HM, Kranges, K_min, sr_F, sr_expt, sr_ζ, q_N1, sr_Icout, sr_expNout, sr_w = method
    b .= -get_T0(medium)
    Nb = n_boreholes(borefield)

    for i in 1:Nb
        @views b[i] -= dot(sr_w[i], sr_expt[i] .* sr_F[i] .- q_N1[i] .* sr_expNout[i] .* (1 .- sr_expt[i]) ./ sr_ζ[i])
    end

    for target in 1:Nb
        for source in 1:Nb
            isempty(Kranges) && continue
            @inbounds range = Kranges[K_min[source, target]]
            @inbounds @views b[target] -= dot(F[range, source], HM[range, target, source])
        end
    end
end

function reset!(method::NonHistoryMethod) 
    @unpack F, sr_F = method
    F .= 0
    for f in sr_F
        f .= 0
    end
end
