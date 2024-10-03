using .FiniteLineSource: initialize_buffer, Buffer

abstract type QuadGKBuffer end

struct NoBoundaryQuadGKBuffer <: QuadGKBuffer
    buffer::Buffer{Float64}
end
NoBoundaryQuadGKBuffer(s, rb::T) where {T <: Number} = NoBoundaryQuadGKBuffer(Buffer(initialize_buffer(s, rb)))

struct DirichletQuadGKBuffer <: QuadGKBuffer
    buffer::Buffer{Float64}
    image_buffer::Buffer{Float64}
end
DirichletQuadGKBuffer(s, rb::T) where {T <: Number} = DirichletQuadGKBuffer(Buffer(initialize_buffer(s, rb)), Buffer(initialize_buffer(image(s), rb)))

get_buffers(::NoBoundary) = Vector{NoBoundaryQuadGKBuffer}()
get_buffers(::DirichletBoundaryCondition) = Vector{DirichletQuadGKBuffer}()

function add_buffer!(buffers, ::NoBoundary, s, rb)
    push!(buffers, NoBoundaryQuadGKBuffer(s, rb))
end

function add_buffer!(buffers, ::DirichletBoundaryCondition, s, rb)
    push!(buffers, DirichletQuadGKBuffer(s, rb))
end