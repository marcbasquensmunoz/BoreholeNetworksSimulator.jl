
mutable struct ConvolutionMethod{T} <: TimeSuperpositionMethod 
    g::Array{T, 3}
    q::Array{T, 2}
end
ConvolutionMethod() = ConvolutionMethod(zeros(0,0,0), zeros(0,0))

function precompute_auxiliaries!(model::ConvolutionMethod; options::SimulationOptions)
    @unpack Nb, Nt, t, borefield, boundary_condition = options
    model.g = zeros(Nb, Nb, Nt)
    model.q = zeros(Nb, Nt)
    compute_response!(model.g, borefield.medium, borefield, boundary_condition, t)
    return model
end

function update_auxiliaries!(method::ConvolutionMethod, X, borefield::Borefield, step)
    Nb = n_boreholes(borefield)
    method.q[:, step] = @view X[3Nb+1:end, step] 
end

function method_coeffs!(M, method::ConvolutionMethod, borefield::Borefield, ::BoundaryCondition)
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