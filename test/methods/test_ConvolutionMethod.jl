import BoreholeNetworksSimulator: method_coeffs!, method_b!, precompute_auxiliaries!, update_auxiliaries!
import BoreholeNetworksSimulator: MediumMock, BorefieldMock, ConstraintMock, BoundaryConditionMock, FluidMock

@testset "test_ConvolutionMethod_auxiliaries" begin
    method = ConvolutionMethod()
    Nb = 2
    Nt = 10
    g = 4.5
    borefield = BorefieldMock(Nb=Nb)
    medium = MediumMock(g=g)
    constraint = ConstraintMock()
    options = SimulationOptions(
        method=method,
        constraint=constraint,
        borefield=borefield,
        medium=medium,
        fluid=FluidMock(),
        Δt=3600*24*30.,
        Nt=Nt,
        configurations=[]
    )
    precompute_auxiliaries!(method, options)

    @test size(method.g) == (Nb, Nb, Nt)
    @test method.g == g .* ones(Nb, Nb, Nt)
    @test size(method.q) == (Nb, Nt)
    @test sum(abs.(method.q)) == 0

    X = zeros(4Nb, Nt)
    for i in 1:Nt
        X[3Nb+1:end, i] .= i 
    end
    for i in 1:Nt
        update_auxiliaries!(method, X, borefield, i)
    end
    @test method.q == X[3Nb+1:end, :]
end

@testset "test_ConvolutionMethod_method_coeffs!" begin
    Nb = 2
    Nt = 10
    method = ConvolutionMethod()
    g = 4.5
    borefield = BorefieldMock(Nb=Nb)
    medium = MediumMock(g=g)
    constraint = ConstraintMock()
    boundary_condition = BoundaryConditionMock()
    options = SimulationOptions(
        method=method,
        constraint=constraint,
        borefield=borefield,
        medium=medium,
        fluid=FluidMock(),
        boundary_condition=boundary_condition,
        Δt=3600*24*30.,
        Nt=Nt,
        configurations=[]
    )
    precompute_auxiliaries!(method, options)

    M = zeros(Nb, 4Nb)
    method_coeffs!(M, method, borefield, medium, boundary_condition)

    expected = [
        (1, 2Nb+1, -1.), (1, 3Nb+1, g), (1, 3Nb+2, g),
        (2, 2Nb+2, -1.), (2, 3Nb+1, g), (2, 3Nb+2, g)
        ]
    @test test_sparse_matrix(M, expected)
end


@testset "test_ConvolutionMethod_method_b!" begin
    Nb = 2
    Nt = 10
    method = ConvolutionMethod()
    g = zeros(Nb, Nb, Nt)

    for i in 1:Nt
        g[:, :, i] .= i
    end

    borefield = BorefieldMock(Nb=Nb)
    medium = MediumMock(g=g, T0=10.)
    constraint = ConstraintMock()
    boundary_condition = BoundaryConditionMock()
    options = SimulationOptions(
        method=method,
        constraint=constraint,
        borefield=borefield,
        medium=medium,
        fluid=FluidMock(),
        boundary_condition=boundary_condition,
        Δt=3600*24*30.,
        Nt=Nt,
        configurations=[]
    )
    precompute_auxiliaries!(method, options)
    method.q = ones(Nb, Nt)

    b = zeros(Nb)

    method_b!(b, method, borefield, medium, 1)
    @test b == -10 .* ones(Nb)

    method_b!(b, method, borefield, medium, 4)
    @test b == -16 .* ones(Nb)

    method_b!(b, method, borefield, medium, 7)
    @test b == -22 .* ones(Nb)

    method_b!(b, method, borefield, medium, 10)
    @test b == -28 .* ones(Nb)
end
