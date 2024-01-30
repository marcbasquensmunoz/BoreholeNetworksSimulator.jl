mutable struct ConvolutionGroundModel{T} <: GroundModel 
    g::Array{T, 3}
    Δq::Array{T, 2}
    T0::T
end
function ConvolutionGroundModel(;T0, parameters)
    @unpack Nb, Nt, Ns = parameters
    ConvolutionGroundModel(zeros(Nb, Nb, Nt), zeros(Nt, Ns), T0)
end

function precompute_auxiliaries!(model::ConvolutionGroundModel, borefield::Borefield, t) 
    compute_response!(model.g, borefield.medium, borefield, t)
end

function update_auxiliaries!(model::ConvolutionGroundModel, X, borefield::Borefield, step)
    Nb = borehole_amount(borefield)
    model.Δq[step, :] = @view X[step, 3Nb+1:end] 
end

function ground_model_coeffs!(M, model::ConvolutionGroundModel, borefield::Borefield)
    Nb = borehole_amount(borefield)
    Ns = segment_amount(borefield)
    M[1:Ns, 3Nb+1:3Nb+Ns] = @view model.g[:,:,1]
    for i in 1:Ns
        bh = where_is_segment(borefield, i)
        M[i, 2Nb + bh] = -1
    end
end

function ground_model_b!(b, model::ConvolutionGroundModel, borefield::Borefield, step)
    Ns = segment_amount(borefield)
    b .= -model.T0 
    for i in 1:Ns
        for j = 1:Ns
            for k = 2:step
                b[i] -= model.Δq[step - k + 1, j] * model.g[j, i, k]
            end
        end
    end
end