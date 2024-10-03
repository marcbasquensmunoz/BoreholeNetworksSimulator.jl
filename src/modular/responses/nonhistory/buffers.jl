using .FiniteLineSource: initialize_buffer

abstract type QuadGKBuffer end

struct NoBoundaryQuadGKBuffer <: QuadGKBuffer
    buffer::Vector{QuadGK.Segment{Float64, Float64, Float64}}
end

struct DirichletQuadGKBuffer <: QuadGKBuffer
    buffer::Vector{QuadGK.Segment{Float64, Float64, Float64}}
    image_buffer::Vector{QuadGK.Segment{Float64, Float64, Float64}}
end

get_buffers(::NoBoundary) = Vector{NoBoundaryQuadGKBuffer}()
get_buffers(::DirichletBoundaryCondition) = Vector{DirichletQuadGKBuffer}()

function add_buffer!(buffers, ::NoBoundary, s, rb)
    push!(buffers, NoBoundaryQuadGKBuffer(initialize_buffer(s, rb)))
end

function add_buffer!(buffers, ::DirichletBoundaryCondition, s, rb)
    buffer = DirichletQuadGKBuffer(initialize_buffer(s, rb), initialize_buffer(image(s), rb))
    push!(buffers, buffer)
end