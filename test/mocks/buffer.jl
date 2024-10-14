using .BoreholeNetworksSimulator: QuadGKBuffer

struct QuadGKBufferMock <: QuadGKBuffer end

BoreholeNetworksSimulator.get_buffers(::BoundaryConditionMock) = Vector{QuadGKBufferMock}()
BoreholeNetworksSimulator.add_buffer!(buffers, ::BoundaryConditionMock, s, rb) = push!(buffers, QuadGKBufferMock())
BoreholeNetworksSimulator.weights(::BoundaryConditionMock, setup::SetupMock, params, dp, containers, ::QuadGKBufferMock; atol, rtol) = setup.weights
