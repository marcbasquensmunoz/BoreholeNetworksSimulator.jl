
get_Nb(containers) = Int(length(containers.b)/4)

get_Tfin(containers)                      = containers.X[get_Tfin_range(containers), :]
get_Tfin(containers, borehole::Int)       = containers.X[2 * borehole - 1, :]
get_Tfin(containers, branch::Vector{Int}) = containers.X[2 .* branch .- 1, :]

get_Tfout(containers)                      = containers.X[get_Tfout_range(containers), :]
get_Tfout(containers, borehole::Int)       = containers.X[2 * borehole, :]
get_Tfout(containers, branch::Vector{Int}) = containers.X[2 .* branch, :]

get_Tb(containers)                      = containers.X[get_Tb_range(containers), :]
get_Tb(containers, borehole::Int)       = containers.X[get_Tb_range_start(containers) + borehole, :]
get_Tb(containers, branch::Vector{Int}) = containers.X[get_Tb_range_start(containers) .+ branch, :]

get_q(containers)                      = containers.X[get_q_range(containers):end,:]
get_q(containers, borehole::Int)       = containers.X[get_q_range_start(containers) + borehole, :]
get_q(containers, branch::Vector{Int}) = containers.X[get_q_range_start(containers) .+ branch, :]

get_mf(containers)                      = containers.mf
get_mf(containers, borehole::Int)       = containers.mf[borehole, :]
get_mf(containers, branch::Vector{Int}) = containers.mf[branch, :]

get_Tfin_range(Nb::Int) = 1:2:2Nb
get_Tfin_range(containers) = get_Tfin_range(get_Nb(containers))
#get_Tfin_range(options) = get_Tfin_range(options.Nb)

get_Tfout_range(Nb::Int) = 2:2:2Nb
get_Tfout_range(containers) = get_Tfout_range(get_Nb(containers))
#get_Tfout_range(options) = get_Tfout_range(options.Nb)

get_Tb_range(Nb::Int) = 2Nb+1:3Nb
get_Tb_range(containers) = get_Tb_range(get_Nb(containers))
#get_Tb_range(options) = get_Tb_range(options.Nb)
get_Tb_range_start(Nb::Int) = 2Nb
get_Tb_range_start(containers) = get_Tb_range_start(get_Nb(containers))

get_q_range(Nb::Int) = 3Nb+1:4Nb
get_q_range(containers) = get_q_range(get_Nb(containers))
#get_q_range(options) = get_q_range(options.Nb)
get_q_range_start(Nb::Int) = 3Nb
get_q_range_start(containers) = get_q_range_start(get_Nb(containers))
