import BoreholeNetworksSimulator: method_coeffs!, method_b!, precompute_auxiliaries!, update_auxiliaries!
import BoreholeNetworksSimulator: MediumMock, BorefieldMock, ConstraintMock, BoundaryConditionMock

@testset "test_NonHistoryMethod_auxiliaries" begin
    n_disc = 20
    method = NonHistoryMethod(n_disc=n_disc)
    Nb = 2
    Nt = 10
    g = 4.5
    α = 1e-6
    λ = 3.
    borefield = BorefieldMock(Nb=Nb, rb=0.1 .* ones(Nb), 
        coordinates=[(0., 0., 0., 100.), (0., 1., 0., 100.)])
    medium = MediumMock(g=g, α=α, λ=λ)
    constraint = ConstraintMock()
    options = SimulationOptions(
        method=method,
        constraint=constraint,
        borefield=borefield,
        medium=medium,
        Δt=3600*24*30.,
        Nt=Nt,
        configurations=[]
    )
    precompute_auxiliaries!(method, options)

    segments_disc = 13 

    @test size(method.F) == ((n_disc+1)*segments_disc, Nb^2)
    @test sum(abs.(method.F)) == 0
    @test length(method.ζ) == (n_disc+1)*segments_disc
    @test size(method.w) == ((n_disc+1)*segments_disc, Nb^2)
    @test length(method.expΔt) == (n_disc+1)*segments_disc

    ζ = load_data("$(@__DIR__)/zeta")
    @test ζ == method.ζ
    w = load_data("$(@__DIR__)/w")
    @test w == method.w
    expΔt = load_data("$(@__DIR__)/expdt")
    @test expΔt == method.expΔt

    X = zeros(4Nb, Nt)
    for i in 1:Nt
        X[3Nb+1:end, i] .= i 
    end
    update_auxiliaries!(method, X, borefield, 1)

    F = load_data("$(@__DIR__)/F")
    @test F == method.F
end

@testset "test_NonHistoryMethod_method_coeffs!" begin
    method = NonHistoryMethod()
    Nb = 2
    Nt = 10
    g = 4.5
    α = 1e-6
    λ = 3.
    W = 1.2
    borefield = BorefieldMock(Nb=Nb, H=100 .* ones(Nb), rb=0.1 .* ones(Nb), 
        coordinates=[(0., 0., 0., 100.), (0., 1., 0., 100.)])
    medium = MediumMock(g=g, α=α, λ=λ, q_coef=W)
    constraint = ConstraintMock()
    boundary_condition = NoBoundary()
    options = SimulationOptions(
        method=method,
        constraint=constraint,
        borefield=borefield,
        medium=medium,
        boundary_condition=boundary_condition,
        Δt=3600*24*30.,
        Nt=Nt,
        configurations=[]
    )
    precompute_auxiliaries!(method, options)

    M = zeros(Nb, 4Nb)
    method_coeffs!(M, method, borefield, medium, boundary_condition)

    expected = [
        (1, 2Nb+1, 1.), (1, 3Nb+1, -W), (1, 3Nb+2, -W),
        (2, 2Nb+2, 1.), (2, 3Nb+1, -W), (2, 3Nb+2, -W)
        ]
    @test test_sparse_matrix(M, expected)
end

@testset "test_NonHistoryMethod_method_b!" begin
    method = NonHistoryMethod()
    Nb = 2
    Nt = 10
    g = 4.5
    α = 1e-6
    λ = 3.
    W = 1.2
    borefield = BorefieldMock(Nb=Nb, H=100 .* ones(Nb), rb=0.1 .* ones(Nb), 
        coordinates=[(0., 0., 0., 100.), (0., 1., 0., 100.)])
    medium = MediumMock(g=g, α=α, λ=λ, q_coef=W)
    constraint = ConstraintMock()
    boundary_condition = NoBoundary()
    options = SimulationOptions(
        method=method,
        constraint=constraint,
        borefield=borefield,
        medium=medium,
        boundary_condition=boundary_condition,
        Δt=3600*24*30.,
        Nt=Nt,
        configurations=[]
    )
    precompute_auxiliaries!(method, options)
    X = zeros(4Nb, Nt)
    for i in 1:Nt
        X[3Nb+1:end, i] .= i 
    end
    update_auxiliaries!(method, X, borefield, 1)

    b = zeros(Nb)
    method_b!(b, method, borefield, medium, 1)

    @test b == 0.03475387541039184 .* ones(Nb)
end
