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

    rpi = 0.02                          # inner pipe radius
    λp = 0.42                           # pipe material conductivity
    rpo = 0.0223                        # outer pipe radius
    dpw = rpo-rpi                       # pipe thickness
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
get_rp(bh::SingleUPipeBorehole{T}) where {T <: Real} = bh.rpi
get_default_hp(bh::SingleUPipeBorehole{T}) where {T <: Real} = bh.hp
get_n_segments(bh::SingleUPipeBorehole) = bh.n_segments


function uniform_Tb_coeffs(borehole::SingleUPipeBorehole, λ, mass_flow, Tref, fluid)
    @unpack λg, λp, rb, rpi, rpo, dpw, H, R_cache, A, method, exp_cache = borehole

    if mass_flow == 0.
        return 0., -1., 1.
    end 

    hp = heat_transfer_coefficient(mass_flow, Tref, borehole, fluid)
    Rhp = 1/(2*π*rpi*hp)

    if iszero(R_cache)
        x1, y1 = borehole.pipe_position[1]
        x2, y2 = borehole.pipe_position[2]

        Rp = 1/(2*π*λp)*log(rpo/rpi)

        k = (λg - λ)/(λg + λ)
        r2_1 = x1^2+y1^2
        r2_2 = x2^2+y2^2
        r2 = (x1 - x2)^2 + (y1 - y2)^2
        d12 = sqrt( (1 - r2_1 / rb^2) * (1 - r2_2 / rb^2) + r2 / rb^2 )

        R11 =  1/(2*π*λg) * (log(rb/rpo) - k * log(1 - r2_1 / rb^2) ) + Rp
        R12 = -1/(2*π*λg) * (log(sqrt(r2) / rb ) + k * log(d12))
        R22 =  1/(2*π*λg) * (log(rb/rpo) - k * log(1 - r2_2 / rb^2) ) + Rp 

        R_cache = @SMatrix [R11 R12; R12 R22]
    end
    R11 = R_cache[1, 1] + Rhp
    R12 = R_cache[1, 2]
    R22 = R_cache[2, 2] + Rhp

    mfc = mass_flow * cpf(fluid)
    det = (R11*R22 - R12^2)

    β1β2 = (R22 - R12 + R11 - R12) / (det * mfc)
    β12 = R12 / (det * mfc)

    γ = sqrt(β1β2^2 / 4 + β1β2 * β12)

    L = β1β2/(2γ) * tanh(γ * H)

    k_i = (1-L)/(1+L)
    return k_i, -1., 1 - k_i
end
