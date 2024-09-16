"""
    SingleUPipeBorehole{T <: Real} <: Borehole @deftype T
    SingleUPipeBorehole(H, D)

Model a borehole with a single U-pipe with burial depth `D` and length `H`.

# Arguments
- `λg = 2.5`: grout conductivity
- `Cg = 2000. * 1550.`: grout capacity
- `αg = λg/Cg`: grout thermal diffusivity
- `rp = 0.02`: pipe radius
- `λp = 0.42`: pipe material conductivity
- `dpw = 0.0023`: pipe thickness
- `rpo = rp + dpw `: equivalent pipe radius
- `hp = 725.`: heat transfer coefficient fluid to pipe
- `pipe_position::NTuple{2, Tuple{T, T}} = [(0.03, 0.0), (-0.03, 0.0)]`: positions of the downward and upward branches of the pipe. (0, 0) represents the center of the borehole.
- `rb = 0.115/2`: borehole radius
"""
@with_kw struct SingleUPipeBorehole{T <: Real} <: Borehole @deftype T
    λg = 2.5                            # grout conductivity
    Cg = 2000. * 1550.                  # grout capacity
    αg = λg/Cg	                        # grout thermal diffusivity

    rp = 0.02                           # pipe radius
    λp = 0.42                           # pipe material conductivity
    dpw = 0.0023                        # pipe thickness
    rpo = rp + dpw                      # equivalent pipe radius
    hp = 725.                           # heat transfer coefficient fluid to pipe ?
    pipe_position::NTuple{2, Tuple{T, T}} = ((0.03, 0.0), (-0.03, 0.0))
        
    rb = 0.115/2                        # borehole radius
    H                                   # length of the borehole
    D                                   # burial depth of the borehole

    n_segments::Int = 1

    R_cache::SMatrix{2, 2, T, 4} = @SMatrix zeros(2, 2)
    A::MMatrix{2, 2, T, 4} = @MMatrix zeros(2, 2)
    method::ExpMethodHigham2005 = ExpMethodHigham2005(false)
    exp_cache::Tuple{Vector{MMatrix{2, 2, T, 4}}, Vector{T}} = ExponentialUtilities.alloc_mem(A, method)
end

get_H(bh::SingleUPipeBorehole{T}) where {T <: Real} = bh.H
get_D(bh::SingleUPipeBorehole{T}) where {T <: Real} = bh.D
get_h(bh::SingleUPipeBorehole{T}) where {T <: Real} = bh.H / bh.n_segments
get_rb(bh::SingleUPipeBorehole{T}) where {T <: Real} = bh.rb
get_rp(bh::SingleUPipeBorehole{T}) where {T <: Real} = bh.rp
get_default_hp(bh::SingleUPipeBorehole{T}) where {T <: Real} = bh.hp
get_n_segments(bh::SingleUPipeBorehole) = bh.n_segments


function uniform_Tb_coeffs(borehole::SingleUPipeBorehole, λ, mass_flow, Tref, fluid)
    @unpack λg, λp, rb, rp, rpo, dpw, H, R_cache, A, method, exp_cache = borehole

    hp = heat_transfer_coefficient(mass_flow, Tref, borehole, fluid.name)

    if iszero(R_cache)
        x1, y1 = borehole.pipe_position[1]
        x2, y2 = borehole.pipe_position[2]

        Rp = 1/(2*π*λp)*log(rp/(rp-dpw))

        d12 = sqrt( (1 - (x1^2+y1^2) / rb^2) * (1 - (x2^2+y2^2) / rb^2) + ( (x1 - x2)^2 + (y1 - y2)^2) / rb^2 )

        R11 =  1/(2*π*λg) * (log(rb/rpo) - (λg - λ)/(λg + λ) * log(1 - (x1^2 + y1^2) / rb^2) ) + Rp
        R12 = -1/(2*π*λg) * (log(( (x1 - x2)^2 + (y1 - y2)^2) / rb^2 ) + (λg - λ)/(λg + λ) * log(d12))
        R22 =  1/(2*π*λg) * (log(rb/rpo) - (λg - λ)/(λg + λ) * log(1 - (x2^2 + y2^2) / rb^2) ) + Rp 
        R_cache = @SMatrix [R11 R12; R12 R22]
    end

    Rhp = 1/(2*π*rp*hp)

    A .= R_cache
    A[1, 1] += Rhp
    A[2, 2] += Rhp

    C = H / (fluid.cpf * mass_flow)
    In = @SMatrix [-1 0; 0 1]
    A .= C .* In * inv(A)

    exponential!(A, method, exp_cache)

    EoutH = A[1, 2] - A[2, 2]
    EinH  = A[2, 1] - A[1, 1]      
    return EinH, -EoutH, EoutH - EinH
end
