
mutable struct ConvolutionMethod{T} <: TimeSuperpositionMethod 
    g::Array{T, 3}
    q::Array{T, 2}
end
function ConvolutionMethod(;parameters, borefield)
    @unpack Nb, Nt, Ns = parameters
    model = ConvolutionMethod(zeros(Nb, Nb, Nt), zeros(Ns, Nt))
    compute_response!(model.g, borefield.medium, borefield, parameters.t)
    return model
end

function update_auxiliaries!(method::ConvolutionMethod, X, borefield::Borefield, step)
    Nb = n_boreholes(borefield)
    @show X[3Nb+1:end, step] 
    method.q[:, step] = @view X[3Nb+1:end, step] 
end

function method_coeffs!(M, method::ConvolutionMethod, borefield::Borefield)
    Nb = n_boreholes(borefield)
    Ns = n_segments(borefield)
    M[1:Ns, 3Nb+1:3Nb+Ns] = @view method.g[:,:,1]
    for i in 1:Ns
        bh = where_is_segment(borefield, i)
        M[i, 2Nb + bh] = -1
    end
end

function method_b!(b, method::ConvolutionMethod, borefield::Borefield, step)
    Ns = n_segments(borefield)
    b .= -get_T0(borefield)

    for k = 1:step-1
        for i in 1:Ns
            @views @inbounds b[i] -= dot(method.q[:, k], method.g[:, i, step - k + 1] - method.g[:, i, step - k])
        end
    end
end