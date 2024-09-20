module BNSPlots

using WGLMakie

include("results.jl")
include("borefield.jl")
include("monitor.jl")

export monitor
export plot_borefield

end
