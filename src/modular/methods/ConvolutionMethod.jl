
mutable struct ConvolutionMethod{T} <: Method 
    g::Array{T, 3}
    Δq::Array{T, 2}
end
function ConvolutionMethod(;parameters, borefield)
    @unpack Nb, Nt, Ns = parameters
    model = ConvolutionMethod(zeros(Nb, Nb, Nt), zeros(Ns, Nt))
    precompute_auxiliaries!(model, borefield, parameters.t)
    return model
end

function precompute_auxiliaries!(method::ConvolutionMethod, borefield::Borefield, t) 
    compute_response!(method.g, borefield.medium, borefield, t)
end

function update_auxiliaries!(method::ConvolutionMethod, X, current_Q, borefield::Borefield, step)
    Nb = borehole_amount(borefield)
    method.Δq[:, step] = @view X[3Nb+1:end, step] 
end

function method_coeffs!(M, method::ConvolutionMethod, borefield::Borefield)
    Nb = borehole_amount(borefield)
    Ns = segment_amount(borefield)
    M[1:Ns, 3Nb+1:3Nb+Ns] = @view method.g[:,:,1]
    for i in 1:Ns
        bh = where_is_segment(borefield, i)
        M[i, 2Nb + bh] = -1
    end
end

function method_b!(b, method::ConvolutionMethod, borefield::Borefield, step, current_Q)
    Ns = segment_amount(borefield)
    b .= -get_T0(borefield)

    for k = 2:step
        for i in 1:Ns
            @views @inbounds b[i] -= dot(method.Δq[:, step - k + 1], method.g[:, i, k])
        end
    end
end