
function topology_coeffs!(M, operation::BoreholeOperation)
    M .= zero(eltype(M))
    j = 1
    for branch in operation.network.branches
        for i in eachindex(@view branch[1:end-1])
            in = branch[i+1] 
            out = branch[i] 
            M[j, 2*in-1] = 1.
            M[j, 2*out] = -1.
            j += 1
        end
    end
end

function topology_b!(b, ::BoreholeOperation) end