mutable struct ConvolutionMethod{T} <: Method 
    g::Array{T, 3}
    Δq::Array{T, 2}
    T0::T
end
function ConvolutionMethod(;T0, parameters, borefield)
    @unpack Nb, Nt, Ns = parameters
    model = ConvolutionMethod(zeros(Nb, Nb, Nt), zeros(Nt, Ns), T0)
    precompute_auxiliaries!(model, borefield, parameters.t)
    return model
end

function precompute_auxiliaries!(method::ConvolutionMethod, borefield::Borefield, t) 
    compute_response!(method.g, borefield.medium, borefield, t)
end

function update_auxiliaries!(method::ConvolutionMethod, X, borefield::Borefield, step)
    Nb = borehole_amount(borefield)
    method.Δq[step, :] = @view X[step, 3Nb+1:end] 
end

function method_coeffs!(M, method::ConvolutionMethod, borefield::Borefield)
    Nb = borehole_amount(borefield)
    Ns = segment_amount(borefield)
    M[1:Ns, 3Nb+1:3Nb+Ns] = method.g[:,:,1]
    for i in 1:Ns
        bh = where_is_segment(borefield, i)
        M[i, 2Nb + bh] = -1
    end
end

function method_b!(b, method::ConvolutionMethod, borefield::Borefield, step)
    Ns = segment_amount(borefield)
    b .= -method.T0 
    for i in 1:Ns
        for j = 1:Ns
            for k = 2:step
                b[i] -= method.Δq[step - k + 1, j] * method.g[j, i, k]
            end
        end
    end
end