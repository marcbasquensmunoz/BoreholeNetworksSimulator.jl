module BNSPlots

using WGLMakie
using GraphMakie
using Graphs

include("results.jl")
include("borefield.jl")
include("monitor.jl")

export monitor
export plot_borefield

end
