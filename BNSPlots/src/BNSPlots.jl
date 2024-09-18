module BNSPlots

using WGLMakie

include("results.jl")
include("borefield.jl")
include("monitor_branch.jl")

export monitor_branch
export plot_borefield

end
