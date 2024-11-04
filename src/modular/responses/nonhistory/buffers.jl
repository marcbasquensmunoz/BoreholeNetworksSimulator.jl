using .FiniteLineSource: initialize_buffer, Buffer

abstract type QuadGKBuffer end

struct SimpleQuadGKBuffer <: QuadGKBuffer
    buffer::Buffer{Float64}
end
SimpleQuadGKBuffer(s, rb::T) where {T <: Number} = SimpleQuadGKBuffer(Buffer(initialize_buffer(s, rb)))

struct ImageQuadGKBuffer <: QuadGKBuffer
    buffer::Buffer{Float64}
    image_buffer::Buffer{Float64}
end
ImageQuadGKBuffer(s, rb::T) where {T <: Number} = ImageQuadGKBuffer(Buffer(initialize_buffer(s, rb)), Buffer(initialize_buffer(image(s), rb)))

get_buffers(::NoBoundary) = Vector{SimpleQuadGKBuffer}()
get_buffers(::DirichletBoundaryCondition) = Vector{ImageQuadGKBuffer}()
get_buffers(::NeumannBoundaryCondition) = Vector{ImageQuadGKBuffer}()

function add_buffer!(buffers, ::NoBoundary, s, rb)
    push!(buffers, SimpleQuadGKBuffer(s, rb))
end

function add_buffer!(buffers, ::DirichletBoundaryCondition, s, rb)
    push!(buffers, ImageQuadGKBuffer(s, rb))
end

function add_buffer!(buffers, ::NeumannBoundaryCondition, s, rb)
    push!(buffers, ImageQuadGKBuffer(s, rb))
end