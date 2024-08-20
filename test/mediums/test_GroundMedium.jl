import BoreholeNetworksSimulator: GroundMedium, get_λ, get_α, get_T0, compute_response!

@testset "test_GroundMedium" begin
    λ = rand()
    α = rand()
    T0 = rand()
    medium = GroundMedium(λ=λ, α=α, T0=T0)

    @test get_λ(medium) == λ
    @test get_α(medium) == α
    @test get_T0(medium) == T0
end

@testset "test_GroundMedium_response_NoBoundary" begin
    λ = 3.
    α = 1e-6
    T0 = 10.
    medium = GroundMedium(λ=λ, α=α, T0=T0)

    r = 1.
    t = 3600*30.

    borehole = SingleUPipeBorehole(H=10., D=0.)
    borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=[(0., 0.), (0., r)])
    boundary_condition = NoBoundary()

    g = zeros(2, 2)
    compute_response!(g, medium, borefield, boundary_condition, t) 

    expected = [0.11246425360481815 0.0008175747207000824; 0.0008175747207000824 0.11246425360481815]
    @test isapprox(expected, g)
end

@testset "test_GroundMedium_response_Dirichlet" begin
    λ = 3.
    α = 1e-6
    T0 = 10.
    medium = GroundMedium(λ=λ, α=α, T0=T0)

    r = 1.
    t = 3600*30.

    borehole = SingleUPipeBorehole(H=10., D=0.)
    borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=[(0., 0.), (0., r)])
    boundary_condition = DirichletBoundaryCondition()

    g = zeros(2, 2)
    compute_response!(g, medium, borefield, boundary_condition, t) 

    expected = [0.1116256193603268 0.0008037621082085207; 0.0008037621082085207 0.1116256193603268]
    @test isapprox(expected, g)
end
