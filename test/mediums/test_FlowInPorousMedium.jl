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
