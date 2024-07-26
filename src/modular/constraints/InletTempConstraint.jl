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

function constant_InletTempConstraint(T_in::Vector{N}, Nt) where {N <: Number}
    Nbr = length(T_in)
    T = zeros(Nbr, Nt)
    for i in 1:Nt
        T[:, i] .= T_in
    end
    HeatLoadConstraint(T)
end

function constraints_coeffs!(M, ::InletTempConstraint, operation)
    M .= 0.
    for (i, borehole) in enumerate(first_boreholes(operation.network))
        M[i, 2*borehole - 1] = 1.
    end
end

function constraints_b!(b, constraint::InletTempConstraint, operation, step)
    b .= constraint.T_in[:, step]
end
