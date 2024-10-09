
function weights(::NoBoundary, setup, params::Constants, dp, containers, buffer::NoBoundaryQuadGKBuffer)
    precompute_coefficients(setup, params=params, dp=dp, containers=containers, buffer=buffer.buffer)
end

function weights(::DirichletBoundaryCondition, setup, params::Constants, dp, containers, buffer::DirichletQuadGKBuffer)
    image_setup = image(setup)
    w1 = precompute_coefficients(setup, params=params, dp=dp, containers=containers, buffer=buffer.buffer)
    w2 = precompute_coefficients(image_setup, params=params, dp=dp, containers=containers, buffer=buffer.image_buffer)
    w1 - w2
end
