"""
    ConvolutionMethod{T} <: TimeSuperpositionMethod 
    ConvolutionMethod()

Use the naïve convolution to compute the thermal response between boreholes. 
It should be initialized without arguments, but it contains the variables:
- `g` stores the unit response between each pair of boreholes evaluated at each time of the simulation. 
It should be precomputed with [`initialize`](@ref).
- `q` stores the heat extraction in each borehole at each time step. It is filled as the simulation runs. 
"""
mutable struct ConvolutionMethod{T} <: TimeSuperpositionMethod 
    g::Array{T, 3}
    q::Array{T, 2}
    aux::Vector{T}
end
ConvolutionMethod() = ConvolutionMethod(zeros(0,0,0), zeros(0,0), zeros(0))

function precompute_auxiliaries!(method::ConvolutionMethod, options)
    @unpack Nb, Nt, t, borefield, medium, boundary_condition = options
    method.g = zeros(Nb, Nb, Nt)
    method.q = zeros(Nb, Nt)
    method.aux = zeros(Nb)
    compute_response!(method.g, medium, borefield, boundary_condition, t)
    return method
end

function update_auxiliaries!(method::ConvolutionMethod, X, borefield, step)
    Nb = n_boreholes(borefield)
    method.q[:, step] = @view X[3Nb+1:end, step] 
end

function method_coeffs!(M, method::ConvolutionMethod, borefield, medium, boundary_condition)
    Nb = n_boreholes(borefield)
    Ns = n_segments(borefield)
    M[1:Ns, 3Nb+1:3Nb+Ns] = @view method.g[:,:,1]
    for i in 1:Ns
        bh = where_is_segment(borefield, i)
        M[i, 2Nb + bh] = -1
    end
end

function method_b!(b, method::ConvolutionMethod, borefield, medium, step)
    @unpack g, q, aux = method

    Ns = n_segments(borefield)
    b .= -get_T0(medium)

    for k in 1:step-1
        for i in 1:Ns
            @views @. aux = g[:, i, step - k + 1] - g[:, i, step - k]
            @views @inbounds b[i] -= dot(q[:, k], aux)
        end
    end
end