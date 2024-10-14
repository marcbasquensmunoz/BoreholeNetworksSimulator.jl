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
