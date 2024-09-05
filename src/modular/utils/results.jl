
Nb(containers::SimulationContainers) = Int(length(containers.b)/4)

Tfin(containers::SimulationContainers)                      = containers.X[Tfin_range(containers), :]
Tfin(containers::SimulationContainers, borehole::Int)       = containers.X[2 * borehole - 1, :]
Tfin(containers::SimulationContainers, branch::Vector{Int}) = containers.X[2 .* branch .- 1, :]

Tfout(containers::SimulationContainers)                      = containers.X[Tfout_range(containers), :]
Tfout(containers::SimulationContainers, borehole::Int)       = containers.X[2 * borehole, :]
Tfout(containers::SimulationContainers, branch::Vector{Int}) = containers.X[2 .* branch, :]

Tb(containers::SimulationContainers)                      = containers.X[Tb_range(containers), :]
Tb(containers::SimulationContainers, borehole::Int)       = containers.X[Tb_range_start(containers) + borehole, :]
Tb(containers::SimulationContainers, branch::Vector{Int}) = containers.X[Tb_range_start(containers) .+ branch, :]

q(containers::SimulationContainers)                      = containers.X[q_range(containers):end,:]
q(containers::SimulationContainers, borehole::Int)       = containers.X[q_range_start(containers) + borehole, :]
q(containers::SimulationContainers, branch::Vector{Int}) = containers.X[q_range_start(containers) .+ branch, :]

Tfin_range(Nb) = 1:2:2Nb
Tfin_range(containers::SimulationContainers) = Tfin_range(Nb(containers))
Tfin_range(options::SimulationOptions) = Tfin_range(options.Nb)

Tfout_range(Nb) = 2:2:2Nb
Tfout_range(containers::SimulationContainers) = Tfout_range(Nb(containers))
Tfout_range(options::SimulationOptions) = Tfout_range(options.Nb)

Tb_range(Nb) = 2Nb+1:3Nb
Tb_range(containers::SimulationContainers) = Tb_range(Nb(containers))
Tb_range(options::SimulationOptions) = Tb_range(options.Nb)
Tb_range_start(Nb) = 2Nb
Tb_range_start(containers::SimulationContainers) = Tb_range_start(Nb(containers))

q_range(Nb) = 3Nb+1:4Nb
q_range(containers::SimulationContainers) = q_range(Nb(containers))
q_range(options::SimulationOptions) = q_range(options.Nb)
q_range_start(Nb) = 3Nb
q_range_start(containers::SimulationContainers) = q_range_start(Nb(containers))
