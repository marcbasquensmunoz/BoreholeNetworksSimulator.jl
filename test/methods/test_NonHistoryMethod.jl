import BoreholeNetworksSimulator: method_coeffs!, method_b!, precompute_auxiliaries!, update_auxiliaries!


global const atol = eps()
@testset "test_OriginalNonHistoryMethod_auxiliaries" begin
    n_disc = 20
    segments_disc = 13 

    method = OriginalNonHistoryMethod(n_disc=n_disc)
    Nb = 2
    Nt = 10
    g = 4.5
    α = 1e-6
    λ = 3.

    weights = rand(n_disc+1)

    borefield = BorefieldMock(Nb=Nb, rb=0.1 .* ones(Nb), 
        coordinates=[(0., 0.), (0., 1.)], D = zeros(Nb), H = 100 .* ones(Nb))
    medium = MediumMock(step_response=g, α=α, λ=λ, q_coef=weights)
    constraint = ConstraintMock()
    fluid = FluidMock()
    boundary_condition = BoundaryConditionMock()
    options = SimulationOptions(
        method=method,
        constraint=constraint,
        borefield=borefield,
        medium=medium,
        boundary_condition=boundary_condition,
        fluid=fluid,
        Δt=3600*24*30.,
        Nt=Nt,
        configurations=[]
    )
    precompute_auxiliaries!(method, options)


    @test size(method.F) == ((n_disc+1)*segments_disc, Nb^2)
    @test sum(abs.(method.F)) == 0
    @test length(method.ζ) == (n_disc+1)*segments_disc
    @test size(method.w) == ((n_disc+1)*segments_disc, Nb^2)
    @test length(method.expΔt) == (n_disc+1)*segments_disc

    ζ = load_data("$(@__DIR__)/zeta")
    @test ζ ≈ method.ζ atol=atol

    w = zeros(size(method.w))
    for j in 1:size(method.w)[2]
        for i in 1:segments_disc
            w[(i-1)*(n_disc+1)+1:i*(n_disc+1), j] .= weights
        end
    end

    expΔt = load_data("$(@__DIR__)/expdt")
    @test expΔt ≈ method.expΔt atol=atol

    X = zeros(4Nb, Nt)
    for i in 1:Nt
        X[3Nb+1:end, i] .= i 
    end
    update_auxiliaries!(method, X, borefield, 1)

    F = load_data("$(@__DIR__)/F")
    @test F ≈ method.F atol=atol
end

@testset "test_OriginalNonHistoryMethod_method_coeffs!" begin
    method = OriginalNonHistoryMethod()
    Nb = 2
    Nt = 10
    g = 4.5
    α = 1e-6
    λ = 3.
    W = 1.2
    borefield = BorefieldMock(Nb=Nb, H=100 .* ones(Nb), D = zeros(Nb), rb=0.1 .* ones(Nb), 
        coordinates=[(0., 0.), (0., 1.)])
    medium = MediumMock(step_response=g, α=α, λ=λ, q_coef=W)
    constraint = ConstraintMock()
    boundary_condition = BoundaryConditionMock()
    fluid = FluidMock()
    options = SimulationOptions(
        method=method,
        constraint=constraint,
        borefield=borefield,
        medium=medium,
        fluid=fluid,
        boundary_condition=boundary_condition,
        Δt=3600*24*30.,
        Nt=Nt,
        configurations=[]
    )
    precompute_auxiliaries!(method, options)

    M = zeros(Nb, 4Nb)
    method_coeffs!(M, method, options)

    expected = [
        (1, 2Nb+1, -1.), (1, 3Nb+1, W), (1, 3Nb+2, W),
        (2, 2Nb+2, -1.), (2, 3Nb+1, W), (2, 3Nb+2, W)
        ]
    @test test_sparse_matrix(M, expected)
end

@testset "test_OriginalNonHistoryMethod_method_b!" begin
    method = OriginalNonHistoryMethod()
    Nb = 2
    Nt = 10
    g = 4.5
    α = 1e-6
    λ = 3.
    W = 1.2
    borefield = BorefieldMock(Nb=Nb, H=100 .* ones(Nb), D=zeros(Nb), rb=0.1 .* ones(Nb), 
        coordinates=[(0., 0.), (0., 1.)])
    medium = MediumMock(step_response=g, α=α, λ=λ, q_coef=W)
    constraint = ConstraintMock()
    boundary_condition = NoBoundary()
    fluid = FluidMock()
    options = SimulationOptions(
        method=method,
        constraint=constraint,
        borefield=borefield,
        medium=medium,
        fluid=fluid,
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

    @test b ≈ -856.5309617326807 .* ones(Nb) atol=atol
end
