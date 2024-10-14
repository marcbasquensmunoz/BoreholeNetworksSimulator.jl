using FiniteLineSource
using QuadGK

struct SetupMock <: FiniteLineSource.Setup 
    weights
end
struct ComputationContainerMock <: FiniteLineSource.ComputationContainers end

BoreholeNetworksSimulator.setup(::Approximation, m::MediumMock, ::Borefield, ::Int64, ::Int64) = SetupMock(m.q_coef)
BoreholeNetworksSimulator.image(s::SetupMock) = SetupMock(s.weights)
BoreholeNetworksSimulator.setup_type(::Approximation, ::MediumMock) = SetupMock
FiniteLineSource.initialize_buffer(::SetupMock, ::Float64) = Vector{QuadGK.Segment{Float64, Float64, Float64}}()
function FiniteLineSource.initialize_containers(::SetupMock, dps)
    N = map(dp -> dp.n, dps)
    unique_n = unique(N)
    map_n = [findfirst(x -> x==n, unique_n) for n in N]
    (map_n, map(n -> ComputationContainerMock(), unique_n))
end
FiniteLineSource.precompute_coefficients(s::SetupMock; params, dp, containers, buffer, atol, rtol) = s.weights
