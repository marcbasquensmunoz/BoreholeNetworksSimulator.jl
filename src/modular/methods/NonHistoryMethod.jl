mutable struct NonHistoryMethod{T} <: Method 
    F::Matrix{T}
    ζ::Vector{T}
    w::Matrix{T}
    expΔt::Vector{T}
end
check_compatibility(borefield::Borefield, ::Constraint, ::NonHistoryMethod) = borefield.medium isa GroundWaterMedium ? NotCompatible("The non-history method is not implemented with ground water flow yet") : Compatible() 

@with_kw struct SegmentToSegment{T <: Number} @deftype T
    D1
    H1
    D2
    H2
    σ

    rmin = σ
    rLR = sqrt(σ^2 + (D2 - D1 - H1)^2     )
    rLL = sqrt(σ^2 + (D2 - D1)^2          )
    rUL = sqrt(σ^2 + (D1 - D2 - H2)^2     )
    rUR = sqrt(σ^2 + (D2 + H2 - D1 - H1)^2)

    r1 = D2 > D1 + H1 ? rLR : rmin
    r2 = H2 > H1 ? (D2 > D1 ? rLL : rmin) : (D2 + H2 > D1 + H1 ? rUR : rmin)
    r3 = H2 > H1 ? (D2 + H2 > D1 + H1 ? rUR : rmin) : (D2 > D1 ? rLL : rmin)
    r4 = D2 + H2 > D1 ? rUL : rmin
end

function NonHistoryMethod(;parameters, borefield)
    @unpack Nb, Nt, Ns, tstep = parameters
    Δt = tstep
    α = get_α(borefield.medium)
    rb = get_rb(borefield, 1) 
    Δt̃ = α*Δt/rb^2

    b = 10.
    n_segment = 20
    f_guess(z) = exp(-10*z^2*Δt) * (1 - exp(-z^2*Δt)) / z
    segments = adaptive_gk_segments(f_guess, 0., b)
    dps = @views [discretization_parameters(s.a, s.b, n_segment) for s in segments]
    ζ  = reduce(vcat, (dp.x for dp in dps)) 
    expΔt = @. exp(-ζ^2 * Δt̃)

    n = length(ζ)
    kg = get_λ(borefield.medium)

    w = zeros(n, Ns*Ns)

    for i in eachindex(Ns)
        for j in eachindex(Ns)
            rb = get_rb(borefield, (i-1)*Ns+j) # Get the right rb
            sts = get_sts(borefield, i, j)
            w[:, (i-1)*Ns+j] =  reduce(vcat, [precompute_coefficients(sts, dp=dp, kg=kg, rb=rb) for dp in dps])
        end
    end

    return NonHistoryMethod(zeros(n, Ns), ζ, w, expΔt)
end

function get_sts(borefield::Borefield, i, j)
    xi, yi, Di, Hi = segment_coordinates(borefield, i)
    xj, yj, Dj, Hj = segment_coordinates(borefield, j)
    σ = i == j ? get_rb(borefield, i) : sqrt((xi-xj)^2 + (yi-yj)^2)
    SegmentToSegment(D1=Di, H1=Hi, D2=Dj, H2=Hj, σ=σ)
end

function update_auxiliaries!(method::NonHistoryMethod, X, current_Q, borefield::Borefield, step)
    @unpack ζ, w, expΔt = method

    for i in eachindex(size(w)[2])
        @. @views w[:, i] = expΔt * w[:, i] + current_Q[step] * (1 - expΔt) / ζ
    end
end

function method_coeffs!(M, method::NonHistoryMethod, borefield::Borefield)
    Nb = borehole_amount(borefield)
    Ns = segment_amount(borefield)
    λ = get_λ(borefield.medium)

    for i in 1:Ns
        for j in 1:Ns
            sts = get_sts(borefield, i, j)
            M[i, 3Nb+j] = q_coef(borefield.medium, method, sts, (i-1)*Ns+j)
        end
    end

    C = 1 / (2*λ^2*π^2)
    for i in 1:Ns
        bh = where_is_segment(borefield, i)
        M[i, 2Nb + bh] = -1/C
    end
end

function method_b!(b, method::NonHistoryMethod, borefield::Borefield, step)
    @unpack w, expΔt, F = method
    Ns = segment_amount(borefield)
    λ = get_λ(borefield.medium)
    C = 1 / (2*λ^2*π^2)
    b .= -get_T0(borefield) / C

    for i in eachindex(b)
        @inbounds b[i] = dot(w[:, i], expΔt .* F[:, i])
    end
end

function q_coef(::GroundMedium, method, sts, i)
    @unpack D1, H1, D2, H2, σ = sts
    β(d) = sqrt(σ^2 + d^2) + d*log(sqrt(σ^2 + d^2) - d)
    constant_integral = π/(2*H2) * (β(D1+H1-D2-H2) + β(D1-D2) - β(D1+H1-D2) - β(D1-D2-H2))
    constant_integral + dot(method.w[:, i], method.ζ)
end

@with_kw struct IntegrationSegment{T <: Number} @deftype T
    a
    b
    I
    E
end
function adaptive_gk_segments(f, a::T, b::T; rtol = sqrt(eps())) where T <: Number
    n_pre = Int(floor((b-a) / 2))
    n = max(n_pre%2 == 0 ? n_pre : n_pre+1, 8)
    x, w, gw = QuadGK.kronrod(n)
    heap = MutableBinaryMaxHeap{IntegrationSegment{T}}()

    I, E = evalrule(f, a, b, x, w, gw)
    push!(heap, IntegrationSegment(a, b, I, E))

    while E > rtol * abs(I)
        s = pop!(heap)
        mid = (s.a+s.b) * convert(eltype(x), 0.5)
        I1, E1 = evalrule(f, s.a, mid, x, w, gw)
        I2, E2 = evalrule(f, mid, s.b, x, w, gw)
        I = I - s.I + I1 + I2
        E = E - s.E + E1 + E2
        push!(heap, IntegrationSegment(s.a, mid, I1, E1))
        push!(heap, IntegrationSegment(mid, s.b, I2, E2))
    end

    extract_all!(heap)
end

function evalrule(f, a, b, x, w, gw)
    s = (b-a) * 0.5

    Ik = f(a + s) * w[end]
    Ig = zero(Ik)

    for i in eachindex(gw)
        fg = f(a + (1+x[2i])*s) + f(a + (1-x[2i])*s)
        fk = f(a + (1+x[2i-1])*s) + f(a + (1-x[2i-1])*s)
        Ig += fg * gw[i]
        Ik += fg * w[2i] + fk * w[2i-1]
    end
    Ik_s, Ig_s = Ik * s, Ig * s
    E = abs(Ik_s - Ig_s)
    return Ik_s, E
end

function discretization_parameters(a, b, n)
    xt, w = gausslegendre(n+1)    
    m = (b-a)/2
    c = (b+a)/2
    x = @. m * xt + c 
    return (x=x, m=m, c=c, n=n, xt=xt, w=w)
end

function precompute_coefficients(sts::SegmentToSegment; rb, kg, dp)
    @unpack m, c, n, xt, w = dp

    C = sqrt(m*π/2) / (2 * π^2 * rb * kg)

    P = zeros(n+1, n+1)
    M = zeros(n+1)

    for k in 0:n
        for s in 1:n+1
            P[s, k+1] = (2k+1) * w[s] * Pl(xt[s], k)
        end
    end

    h_mean_sts, r_min, r_max = mean_sts_evaluation(sts)
    guide(r) = h_mean_sts(r*rb) * besselj(1/2, r) * imag(exp(im*r)) / r^(3/2)
    R̃, wz = adaptive_gk(guide, r_min/rb, r_max/rb)

    f(r̃, k) = h_mean_sts(r̃*rb) * rb * besselj(k+1/2, m * r̃) * imag((im)^k * exp(im*c*r̃)) / r̃^(3/2)

    for k in 0:n
        M[k+1] = dot(f.(R̃, k), wz)
    end

    M .= P * M 
    return C .* M
end

transpose(p::SegmentToSegment) = SegmentToSegment(D1=p.D2, H1=p.H2, D2=p.D1, H2=p.H1, σ=p.σ)

function L(r, params::SegmentToSegment)
    @unpack D1, H1, D2, H2, σ, r1, r2, r3, r4 = params

    if r <= r1
        0.
    elseif r < r2
        r * (D1 - D2 + H1) / sqrt(r^2-σ^2) + r
    elseif r < r3
        r * min(H1, H2) / sqrt(r^2-σ^2)
    elseif r < r4
        r * (D2 - D1 + H2) / sqrt(r^2-σ^2) - r
    else
        0.
    end
end

function mean_sts_evaluation(params::SegmentToSegment)
    paramsT = transpose(params)
    f(r) = L(r, params) + L(r, paramsT)
    h_mean_sts(r) = f(r) / params.H2
    r_min = Base.max(params.r1, paramsT.r1)
    r_max = Base.max(params.r4, paramsT.r4)
    return h_mean_sts, r_min, r_max
end

function adaptive_gk(f, a, b; rtol = sqrt(eps()))
    x, w, gw = QuadGK.kronrod(8)
    heap = MutableBinaryMaxHeap{IntegrationSegment{Float64}}()

    I, E = evalrule(f, a, b, x, w, gw)
    push!(heap, IntegrationSegment(a, b, I, E))

    while E > rtol * abs(I)
        s = pop!(heap)
        mid = (s.a+s.b) * 0.5
        I1, E1 = evalrule(f, s.a, mid, x, w, gw)
        I2, E2 = evalrule(f, mid, s.b, x, w, gw)
        I = I - s.I + I1 + I2
        E = E - s.E + E1 + E2
        push!(heap, IntegrationSegment(s.a, mid, I1, E1))
        push!(heap, IntegrationSegment(mid, s.b, I2, E2))
    end
    extract_nodes_weights(heap, x, w)
end
Base.isless(x::IntegrationSegment{T}, y::IntegrationSegment{T}) where T <: Number = x.E < y.E

function extract_nodes_weights(heap, x, w)
    res_x = zeros(0)
    res_w = zeros(0)

    while !isempty(heap)
        s = pop!(heap)
        xs, ws = rescale_x_w(s.a, s.b, x, w)
        append!(res_x, xs)
        append!(res_w, ws)
    end

    res_x, res_w
end

function rescale_x_w(a, b, x, w)
    m = (b-a) * 0.5
    c = (b+a) * 0.5
    m .* [x; -x[1:end-1]] .+ c, m .* [w; w[1:end-1]]
end