
function weights(::NoBoundary, setup, params::Constants, dp, containers, buffer::SimpleQuadGKBuffer; atol, rtol)
    precompute_coefficients(setup, params=params, dp=dp, containers=containers, buffer=buffer.buffer, atol=atol, rtol=rtol)
end

function weights(::DirichletBoundaryCondition, setup, params::Constants, dp, containers, buffer::ImageQuadGKBuffer; atol, rtol)
    image_setup = image(setup)
    w1 = precompute_coefficients(setup, params=params, dp=dp, containers=containers, buffer=buffer.buffer, atol=atol, rtol=rtol)
    w2 = precompute_coefficients(image_setup, params=params, dp=dp, containers=containers, buffer=buffer.image_buffer, atol=atol, rtol=rtol)
    w1 - w2
end

function weights(::NeumannBoundaryCondition, setup, params::Constants, dp, containers, buffer::ImageQuadGKBuffer; atol, rtol)
    image_setup = image(setup)
    w1 = precompute_coefficients(setup, params=params, dp=dp, containers=containers, buffer=buffer.buffer, atol=atol, rtol=rtol)
    w2 = precompute_coefficients(image_setup, params=params, dp=dp, containers=containers, buffer=buffer.image_buffer, atol=atol, rtol=rtol)
    w1 + w2
end
