import BoreholeNetworksSimulator: GroundMedium, get_λ, get_α, get_T0, compute_response!

@testset "test_FlowInPorousMedium" begin
    λw = rand()
    λs = rand()
    Cw = rand()
    Cs = rand()
    θ = rand()
    Φ = rand()
    ux_in_meterperday = 1e-2
    T0 = rand()
    medium = FlowInPorousMedium(λw=λw, λs=λs, Cw=Cw, Cs=Cs, θ=θ, Φ=Φ, ux_in_meterperday=ux_in_meterperday, T0=T0)

    @test get_λ(medium) == λs * (1-Φ) + λw*Φ
    @test get_α(medium) == (λs * (1-Φ) + λw*Φ) / (Cs * (1-Φ) + Cw*Φ)
    @test get_T0(medium) == T0

    @test medium.ux == ux_in_meterperday/(3600*24)
    @test medium.vt == ux_in_meterperday/(3600*24) * Cw / (Cs * (1-Φ) + Cw*Φ)    
end


@testset "test_FlowInPorousMedium_response_NoBoundary" begin
    λw = 0.6
    λs = 2.
    Cw = 4.18*1e6
    Cs = 1.7*1e6
    θ = 0.
    Φ = 0.2
    ux_in_meterperday = 1e-2
    T0 = 10.
    medium = FlowInPorousMedium(λw=λw, λs=λs, Cw=Cw, Cs=Cs, θ=θ, Φ=Φ, ux_in_meterperday=ux_in_meterperday, T0=T0)

    r = 1.
    t = 3600*30.

    borehole = SingleUPipeBorehole(H=10., D=0.)
    borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=[(0., 0.), (0., r)])
    boundary_condition = NoBoundary()
    approximation = MeanApproximation()

    g = zeros(2, 2)
    compute_response!(g, medium, borefield, boundary_condition, approximation, t) 

    expected = [0.18534392810902836 0.0006201185607419583; 0.0006201185607419583 0.18534392810902836]
    @test isequal(expected, g)
end


@testset "test_FlowInPorousMedium_response_Dirichlet" begin
    λw = 0.6
    λs = 2.
    Cw = 4.18*1e6
    Cs = 1.7*1e6
    θ = 0.
    Φ = 0.2
    ux_in_meterperday = 1e-2
    T0 = 10.
    medium = FlowInPorousMedium(λw=λw, λs=λs, Cw=Cw, Cs=Cs, θ=θ, Φ=Φ, ux_in_meterperday=ux_in_meterperday, T0=T0)

    r = 1.
    t = 3600*30.

    borehole = SingleUPipeBorehole(H=10., D=0.)
    borefield = EqualBoreholesBorefield(borehole_prototype=borehole, positions=[(0., 0.), (0., r)])
    boundary_condition = DirichletBoundaryCondition()
    approximation = MeanApproximation()

    g = zeros(2, 2)
    compute_response!(g, medium, borefield, boundary_condition, approximation, t) 

    expected = [0.1840767782869683 0.0006106982692048766; 0.0006106982692048766 0.1840767782869683]
    @test isequal(expected, g)
end
