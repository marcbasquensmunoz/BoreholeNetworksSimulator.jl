using .BoreholeNetworksSimulator: QuadGKBuffer

struct QuadGKBufferMock <: QuadGKBuffer end

BoreholeNetworksSimulator.get_buffers(::BoundaryConditionMock) = Vector{QuadGKBufferMock}()
BoreholeNetworksSimulator.add_buffer!(buffers, ::BoundaryConditionMock, s, rb) = push!(buffers, QuadGKBufferMock())
BoreholeNetworksSimulator.weights(::BoundaryConditionMock, setup, params, dp, containers, buffer::QuadGKBufferMock) = ones(length(dp.x))

