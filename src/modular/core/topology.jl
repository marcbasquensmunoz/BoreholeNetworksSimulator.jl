
function topology_coeffs!(M, network::BoreholeNetwork, mass_flows)
    M .= zero(eltype(M))
    j = 1
    Nb = n_boreholes(network)

    for bh_in in 1:Nb
        if bh_in == sink(network)
            continue
        end
        parents = inneighbors(network.graph, bh_in)
        filter!(i -> i ≠ source(network), parents)
        for bh_out in parents
            if bh_out == source(network)
                continue
            end
            M[j, 2*bh_in-1] = 1.
            M[j, 2*bh_out] = - mass_flows[bh_out] / mass_flows[bh_in]
        end
        if !isempty(parents)
            j += 1
        end
    end
end

function topology_b!(b, ::BoreholeOperation) end