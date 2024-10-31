"""
    InletTempConstraint(T_in::Matrix{T}){T <: Number} <: Constraint

Constrain the inlet temperature of the first borehole in each branch.

The inlet temperature `T_in` must be a `Matrix`, whose column `i` are the inlet temperatures per branch at the 
time step `i`. The amount of rows of `T_in` must equal to the amount of branches specified 
in [`BoreholeNetwork`](@ref).
"""
struct InletTempConstraint{T <: Number} <: Constraint 
    T_in::Matrix{T}
end

"""
    constant_InletTempConstraint(T_in::Vector{N}, Nt) where {N <: Number}

Convenience initializer for `InletTempConstraint`. It creates a constant inlet temperature constraint through all the `Nt` time steps,
where `T_in` are the inlet temperatures for each branch.
"""
function constant_InletTempConstraint(T_in::Vector{N}, Nt) where {N <: Number}
    Nbr = length(T_in)
    T = zeros(Nbr, Nt)
    for i in 1:Nt
        T[:, i] .= T_in
    end
    InletTempConstraint(T)
end

"""
    uniform_InletTempConstraint(T_in::Vector{N}, Nbr) where {N <: Number}

Convenience initializer for `InletTempConstraint`. It creates a uniform inlet temperature constraint along all branches,
where `T_in` are the inlet temperatures for each time step.
"""
function uniform_InletTempConstraint(T_in::Vector{N}, Nbr) where {N <: Number}
    Nt = length(T_in)
    T = zeros(Nbr, Nt)
    for i in 1:Nbr
        T[i, :] .= T_in
    end
    InletTempConstraint(T)
end

function constraints_coeffs!(M, ::InletTempConstraint, borefield::Borefield, network, mass_flows)
    M .= zero(eltype(M))

    for (i, borehole) in enumerate(first_bhs_in_branch(network))
        M[i, 2*borehole - 1] = one(eltype(M))
        if mass_flows[borehole] == zero(eltype(M))
            M[i, 2*borehole ] = -one(eltype(M))
        end
    end
end

function constraints_b!(b, constraint::InletTempConstraint, network, mass_flows, step)
    for (i, borehole) in enumerate(first_bhs_in_branch(network))
        if mass_flows[borehole] == zero(eltype(M))
            b[i] = zero(eltype(M))
        else
            b[i] = constraint.T_in[i, step]
        end
    end
end
