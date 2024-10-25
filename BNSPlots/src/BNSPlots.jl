module BNSPlots

using WGLMakie
using GraphMakie

include("results.jl")
include("borefield.jl")
include("monitor.jl")

export monitor
export plot_borefield

end
