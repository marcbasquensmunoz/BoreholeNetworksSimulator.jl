
mutable struct ConvolutionGroundModel <: GroundModel 
    g
    Δq
    T0
end
ConvolutionGroundModel(;T0) = ConvolutionGroundModel([], [], T0)

function precompute_auxiliaries!(model::ConvolutionGroundModel, borefield::Borefield, t) 
    coord_source, coord_eval = segment_coordinates(borefield)
    model.g = response(borefield.medium, borefield, coord_source, coord_eval, t)
    model.Δq = zeros(length(t), segment_amount(borefield))
end

function update_auxiliaries!(model::ConvolutionGroundModel, X, borefield::Borefield, step)
    Nb = borehole_amount(borefield)
    model.Δq[step, :] = X[step, 3Nb+1:end] 
end

function ground_model_coeffs!(M, model::ConvolutionGroundModel, borefield::Borefield)
    Nb = borehole_amount(borefield)
    Ns = segment_amount(borefield)
    M[1:Ns, 3Nb+1:3Nb+Ns] = model.g[:,:,1]
    map = segment_map(borefield)
    for i in 1:Ns
        M[i, 2Nb + map[i]] = -1
    end
end

function ground_model_b!(b, model::ConvolutionGroundModel, borefield::Borefield, step)
    Ns = segment_amount(borefield)
    for i in 1:Ns
        b[i] = -model.T0 
        for j = 1:Ns
            for k = 2:step
                b[i] += -model.Δq[step - k + 1, j] * model.g[j, i, k]
            end
        end
    end
end