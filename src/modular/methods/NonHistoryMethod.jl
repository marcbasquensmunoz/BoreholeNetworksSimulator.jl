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
    g::Vector{Vector{T}}
    first_block_buffer::Vector{CircularBuffer{T}}
end
NonHistoryMethod() = NonHistoryMethod(zeros(0), zeros(0, 0), zeros(0), zeros(0), zeros(0), zeros(0, 0, 0), UnitRange{Int}[], UnitRange{Int}[], zeros(Int, 0, 0), zeros(0), zeros(Int, 0), Vector{Float64}[], CircularBuffer{Float64}[])

get_setup(::MidPointApproximation) = SegmentToPoint(H=0., D=0., z=0., σ=0.)
get_setup(::MeanApproximation) = SegmentToSegment(H1=0., D1=0., H2=0., D2=0., σ=0.)

function precompute_auxiliaries!(method::NonHistoryMethod, options)
    constants = Constants(Δt=options.Δt, α=get_α(options.medium), rb=get_rb(options.borefield, 1), kg=get_λ(options.medium))
    containers = AsymptoticContainers(10)
    ϵ = options.atol != 0 ? options.atol : 1e-4

    Nb = n_boreholes(options.borefield)
    sources = [LineSource(segment_coordinates(options.borefield, i)..., get_rb(options.borefield, i)) for i in 1:Nb]

    rep_setup = get_representative(options.approximation, options.boundary_condition, sources)
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
    method.g = blocks.g
    method.first_block_buffer = blocks.load_buffer
end

function update_auxiliaries!(method::NonHistoryMethod, X, borefield, step)
    @unpack F, N, expt, expNin, expNout, qaux, ranges, first_block_buffer = method

    Nb = n_boreholes(borefield)
    K = length(N) - 1
    @views q = X[3Nb+1:4Nb, 1:step]

    for j in 1:Nb
        push!(first_block_buffer[j], q[j, end])
        qaux .= 0.
        for i in 1:K
            qin = step - N[i] + 1 > 0 ? q[j, step - N[i] + 1] : 0.
            qout = step - N[i+1] + 1 > 0 ? q[j, step - N[i+1] + 1] : 0.
            @views @inbounds @. qaux[ranges[i]] = qin * expNin[ranges[i]] - qout * expNout[ranges[i]]
        end
        @views @. F[:, j] = expt * F[:, j] + qaux
    end
end

image_strength(::NoBoundary) = 0.
image_strength(::DirichletBoundaryCondition) = -1.
image_strength(::NeumannBoundaryCondition) = 1.

get_representative(::MidPointApproximation, bc, sources) = FiniteLineSource.get_representative_ltp(sources, image_strength(bc))[1] 
get_representative(::MeanApproximation, bc, sources) = FiniteLineSource.get_representative_ltl(sources, image_strength(bc))[1] 

function method_coeffs!(M, ::NonHistoryMethod, options)
    @unpack borefield, medium, boundary_condition, approximation, Δt, atol, rtol = options
    @unpack λ, α = medium
    Nb = n_boreholes(borefield)
    Ns = n_segments(borefield)

    for i in 1:Ns
        for j in 1:Ns
            rb = get_rb(borefield, i)
            s = setup(approximation, medium, borefield, i, j)
            M[i, 3Nb+j] = response(boundary_condition, s, Constants(α=α, kg=λ, rb=rb), Δt, atol=atol, rtol=rtol)          
        end
    end

    for i in 1:Ns
        bh = where_is_segment(borefield, i)
        M[i, 2Nb + bh] = -1
    end
end

function method_b!(b, method::NonHistoryMethod, borefield, medium, step)
    @unpack F, HM, Kranges, K_min, g, first_block_buffer = method
    b .= -get_T0(medium)
    Nb = n_boreholes(borefield)

    for i in 1:Nb
        N = length(g[i])
        for (k, qq) in enumerate(first_block_buffer[i])
            if k == 1 continue end
            @inbounds b[i] -= qq * (g[i][N-k+2] - g[i][N-k+1])
        end
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
