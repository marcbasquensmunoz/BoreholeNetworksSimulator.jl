
function test_sparse_matrix(M, expected)
    positions = map(x -> (x[1], x[2]), expected)
    values = map(x -> x[3], expected)
    for j in axes(M, 2) 
        for i in axes(M, 1)
            index = findfirst(x -> x[1]==i && x[2]==j, positions)
            if isnothing(index)
                M[i, j] != 0. && return false
            else 
                M[i, j] != values[index] && return false
            end
        end
    end 
    return true
end