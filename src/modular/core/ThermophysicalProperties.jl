
@with_kw struct ThermophysicalProperties{R <: Number}
    T::Vector{R}
    μ::Vector{R}
    ρ::Vector{R}
    cp::Vector{R}
    k::Vector{R}
end
